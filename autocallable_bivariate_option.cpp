#include "utilities.hpp"


class AutocallableBivariate{
    public:
        mutable NormalRandomGenerator g;
        double payOffAutocallableBivariate(const double& S01, const double& Sref1, double  const & S02, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT, const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&)) const;
        double prixAutocallableBivariate(const double& S01, const double& Sref1, double  const & S02, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT,  const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&), const int& Nmc) const;
        void surfaceAutocallableBivariate(const double& S0max1, size_t const & nbS1, const double& Sref1, double  const & S0max2, size_t const & nbS2, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT,  const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&), const int& Nmc) const;


};


double verifBivariate(const std::vector<double>& S1, const double& Sref1, const std::vector<double>& S2, const double& Sref2, const double& B, const double& r, const int& T, size_t const & nbT,  const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&)) {
    bool flag = true;
    size_t i = 0;
    double payoff;
    while (flag && (i < nbT) ){
        if (std::max(S1[i]/Sref1, S2[i]/Sref2) > B){
            flag = false;
            payoff = exp( -r * t[i] ) * Q_list[i];
        }
        ++i;
    }
    if (flag)
        payoff = exp( -r * T ) * rbmt(S1[nbT-1]/Sref1, S2[nbT-1]/Sref2);
    return payoff;
}


double AutocallableBivariate::payOffAutocallableBivariate(const double& S01, const double& Sref1, double  const & S02, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT, const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&)) const{
    std::vector<double> S1(nbT), S2(nbT);
    double b1, b2;
    for (size_t k = 0; k < nbT; ++k){
        b1 = g();
        b2 = rho * b1 + sqrt(1 - rho*rho) * g();
        S1[k] =  S01 * exp( (r - sigma1*sigma1/2) * t[k] + sigma1 * sqrt(t[k]) * b1 );
        S2[k] =  S02 * exp( (r - sigma2*sigma2/2) * t[k] + sigma2 * sqrt(t[k]) * b2 );
    }
    return verifBivariate(S1, Sref1, S2, Sref2, B, r, T, nbT, t, Q_list, rbmt);
}


double AutocallableBivariate::prixAutocallableBivariate(const double& S01, const double& Sref1, double  const & S02, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT, const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&), const int& Nmc) const{
    double prix = 0;
    for (size_t k = 0; k < Nmc; ++k){
        prix += payOffAutocallableBivariate(S01, Sref1, S02, Sref2, rho, B, r, sigma1, sigma2,T, nbT, t, Q_list, rbmt);
    }
    return prix / Nmc;
}
    

double rbmt2(const double& s1, const double& s2){
    return 100 * std::min(s1, s2);
}


void AutocallableBivariate::surfaceAutocallableBivariate(const double& S0max1, size_t const & nbS1, const double& Sref1, double  const & S0max2, size_t const & nbS2, const double& Sref2, const double& rho, const double& B, const double& r, const double& sigma1, const double& sigma2, const int& T, size_t const & nbT, const std::vector<double>& t, const std::vector<double>& Q_list, double (*rbmt)(const double&, const double&), const int& Nmc) const{
    std::vector<double> S0_list1(nbS1 + 1), S0_list2(nbS2 + 1);
    std::vector<std::vector<double>> prix(nbS1 + 1, std::vector<double>(nbS2 + 1));
    double ds1, ds2;
    ds1 = S0max1/nbS1;
    ds2 = S0max2/nbS2;
    for (int k = 0; k <= nbS1; k++)
        S0_list1[k] = k*ds1;        
    for (int k = 0; k <= nbS2; k++)
        S0_list2[k] = k*ds2;
    for (int n = 0; n <= nbS1; n++){
        for (int i = 0; i <= nbS2; i++)
            prix[n][i] = prixAutocallableBivariate(S0_list1[n], Sref1, S0_list2[i], Sref2, rho, B, r, sigma1, sigma2, T, nbT, t, Q_list, rbmt, Nmc);
    }
    export_data_surface(S0_list1, S0_list2, prix, "surfaceAutocallableBivariate");
}


int main(){

    double S01(350), S02(700), S0max1(700), S0max2(1000), K(350), B(1), T(5), r(0.04), rho(0.5), sigma1(0.3), sigma2(0.4), h(0.1), Sref1(400), Sref2(800);
    int Nmc(10000), nbS1(100), nbS2(100), nbT(5);
    
    std::vector<double> Q_list = {110, 120, 130, 140, 150};
    std::vector<double> t = {1, 2, 3, 4, 5};

    AutocallableBivariate bivar;
    
    bivar.surfaceAutocallableBivariate(S0max1, nbS1, Sref1, S0max2, nbS2, Sref2,rho, B, r, sigma1, sigma2, T, nbT, t, Q_list, rbmt2, Nmc);

    return 0;
}