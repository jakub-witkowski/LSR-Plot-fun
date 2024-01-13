#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <thread>
#include <cstring>
#include <cmath>

/* Include ROOT Framework classes */
#include "TCanvas.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLegend.h"
//#include "TFrame.h"

/* Vector holding each line of the CSV file as string */
std::vector<std::string> raw_data{};

/* Vectors for individual data columns */
std::vector<double> depths{};
std::vector<double> ages{};
std::vector<double> lsr_values{};
std::vector<double> lsr_plot_values{};
std::vector<double> lsr_plot_ages{};
std::vector<double> fit_line{};
std::vector<double> fit_line_x{};
std::vector<double> smoothed_lsr_values{};
std::vector<double> smoothed_lsr_plot_values{};

/* array holding the fitting results */
double chi2{};
int ndf{};
double par[10]{};

bool is_data_sorted(std::vector<double> d)
{
    bool result{true};

    for (int i = 1; i < d.size(); i++)
    {
        if (d[i] < d[i-1])
            result = false;
    }

    return result;
}

double compute_polynomial_expression(int deg, double current_value)
{
    double temporary_result{};
    double total = par[0];
    for (int i = 1; i <= deg; i++)
    {
        temporary_result = par[i] * pow(current_value, i);
        total += temporary_result;
    }
    return total;
}

void age_vs_depth_plot(size_t l, std::vector<double> a, std::vector<double> d, int deg)
{
    auto multi = new TMultiGraph();

    std::string p = "pol " + std::to_string(deg);
    auto f1 = new TF1("f1", p.c_str()); 

    auto g1 = new TGraph(l, &a[0], &d[0]);
    g1->SetName("g1");
    g1->SetTitle("Age vs depth, raw");
    g1->SetMarkerColor(4);
    g1->SetMarkerSize(1.25);
    g1->SetMarkerStyle(20);
    g1->Fit(f1, "N");
    g1->GetXaxis()->SetTicks("-");
    g1->GetYaxis()->SetTicks("-");

    chi2 = f1->GetChisquare();
    ndf = f1->GetNDF();
    for (int i = 0; i <= deg; i++)
    {
        par[i] = f1->GetParameter(i);
    }

    for (int i = 0; i < a.size(); i++)
    {
        fit_line.push_back(compute_polynomial_expression(deg, a[i]));
        fit_line_x.push_back(i);
    }

    auto g2 = new TGraph(fit_line_x.size(), &ages[0], &fit_line[0]);
    g2->SetName("g2");
    g2->SetTitle("Polynomial fit");
    g2->SetLineColor(2);
    g2->SetLineWidth(2);
    g2->GetXaxis()->SetTicks("-");
    g2->GetYaxis()->SetTicks("-");    

    multi->Add(g1, "p");
    multi->Add(g2, "l");
    multi->SetName("AvD");
    multi->SetTitle("Age vs depth plot with polynomial smoothing; Age (Ma);");
    multi->GetXaxis()->CenterTitle();
    multi->GetYaxis()->CenterTitle();
    multi->Draw("A RY");

    auto leg = new TLegend(0.5,0.8,0.9,0.9);
    leg->AddEntry("g1","Age model tiepoints","p");
    leg->AddEntry("g2","Polynomial fit","l");
    leg->Draw();

}

void lsr_plot(size_t l, std::vector<double> a, std::vector<double> d)
{
    auto multi = new TMultiGraph();
    
    auto g1 = new TGraph(l, &a[0], &d[0]);
    g1->SetName("g1");
    g1->SetTitle("LSR, raw");
    g1->SetLineColor(4);
    g1->SetLineWidth(2);
    g1->GetXaxis()->SetTicks("-");
    g1->GetYaxis()->SetTicks("-");

    auto g2 = new TGraph(l, &a[0], &smoothed_lsr_plot_values[0]);
    g2->SetName("g2");
    g2->SetTitle("LSR, smoothed");
    g2->SetLineColor(2);
    g2->SetLineWidth(2);
    g2->GetXaxis()->SetTicks("-");
    g2->GetYaxis()->SetTicks("-");

    multi->Add(g1, "l");
    multi->Add(g2, "l");
    multi->SetName("LSRcomp");
    multi->SetTitle("Raw vs smoothed LSR plot; Age (Ma); Linear sedimentation rate (cm/kyr)");
    multi->GetXaxis()->CenterTitle();
    multi->GetYaxis()->CenterTitle();
    multi->Draw("AL");

    auto leg = new TLegend(0.1,0.8,0.4,0.9);
    leg->AddEntry("g1","LSR, age model","l");
    leg->AddEntry("g2","LSR, smoothed","l");
    leg->Draw();
}

int main(int argc, char** argv)
{
    
    std::string filename{argv[1]};
    std::ifstream input{filename};

    for (std::string line; std::getline(input, line);)
    {
        raw_data.push_back(line);
    }

    /* Separate data on coma, convert substring to double and copy to the respective depths or ages vector */
    for (int i = 0; i < raw_data.size(); i++)
    {
        size_t pos = raw_data[i].find(",");
        depths.push_back(stod(raw_data[i].substr(0, pos)));
        ages.push_back(stod(raw_data[i].substr(pos + 1)));
    }

    /* Test if data is sorted */
    if (!is_data_sorted(depths))
    {
        std::cout << "The data is not properly sorted; inspect the " << filename << " file and try again." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!is_data_sorted(ages))
    {
        std::cout << "The data is not properly sorted; inspect the " << filename << " file and try again." << std::endl;
        exit(EXIT_FAILURE);
    }

    /* calculate LSR values */
    for (int i = 1; i < raw_data.size(); i++)
    {
        lsr_values.push_back(((depths[i]-depths[i-1])*100)/((ages[i]-ages[i-1])*1000));
    }

    /* establish LSR values and age tiepoints for plotting vector */
    for (int i = 0; i < raw_data.size() - 1; i++)
    {
        lsr_plot_values.push_back(lsr_values[i]);
        lsr_plot_values.push_back(lsr_values[i]);
        lsr_plot_ages.push_back(ages[i]);
        lsr_plot_ages.push_back(ages[i+1]);
    }

    auto cnv = new TCanvas("cnv", "cnv", 0, 0, 1200, 800);
    cnv->Divide(2,1);

    cnv->cd(1);
    age_vs_depth_plot(raw_data.size(), ages, depths, atoi(argv[2]));
    
    /* calculate smoothed LSR values */
    for (int i = 1; i < raw_data.size(); i++)
    {
        smoothed_lsr_values.push_back(((fit_line[i]-fit_line[i-1])*100)/((ages[i]-ages[i-1])*1000));
    }

    /* establish LSR values and age tiepoints for plotting vector */
    for (int i = 0; i < raw_data.size() - 1; i++)
    {
        if (lsr_values[i] == 0)
        {
            smoothed_lsr_plot_values.push_back(0);
            smoothed_lsr_plot_values.push_back(0);
        }
        else
        {
            smoothed_lsr_plot_values.push_back(smoothed_lsr_values[i]);
            smoothed_lsr_plot_values.push_back(smoothed_lsr_values[i]);
        }
    }

    cnv->cd(2);
    lsr_plot(lsr_plot_ages.size(), lsr_plot_ages, lsr_plot_values);

    cnv->Update();
    cnv->Modified();

    std::cout << "chi2/ndf: " << chi2/ndf << std::endl;
    cnv->Print("plot.png");

    return 0;
}