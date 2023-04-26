#include "utilities.hpp"


class AutocallableUnivariate{
    public:
        mutable NormalRandomGenerator g;
        double payOffAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT , std::vector<double> const & Q_list, double (*rbmt)(double)) const;
        double prixAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
        void graphAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
        double deltaAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
        void graphDeltaAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
        double gammaAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
        void graphGammaAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const;
};


double verif(std::vector<double> const & S, const double & Sref, const double & B, const int & T, size_t const & nbT, const double & r, std::vector<double> const & t, std::vector<double> const & Q_list, double (*rbmt)(double)) {
    bool flag = true;
    size_t i = 0;
    double payoff;
    while (flag && i < nbT) {
        if (S[i]/Sref > B) {
            flag = false;
            payoff = exp(-r * t[i]) * Q_list[i];
        }
        i++;
    }
    if (flag) {
        payoff = exp(-r * T) * rbmt(S[nbT - 1]/Sref);
    }
    return payoff;
}


double AutocallableUnivariate::payOffAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT , std::vector<double> const & Q_list, double (*rbmt)(double)) const{ 
    std::vector<double> S(nbT), t(nbT);
    for (size_t k = 0; k < nbT; ++k){
        t[k] = (k+1);
        S[k] = S0 * exp( (r - sigma*sigma/2) * t[k] + sigma * sqrt(t[k]) * g());
    }
    return verif(S, Sref, B, T, nbT, r, t, Q_list, rbmt);
}


double AutocallableUnivariate::prixAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    double prix = 0; 
    for (size_t k = 0; k < Nmc; ++k)
        prix += payOffAutocallableUnivariate(S0, Sref, B, r, sigma, T, nbT, Q_list, rbmt);
    return prix / Nmc;
}


double rbmt1(const double s){
    return 100 * s;
}


void AutocallableUnivariate::graphAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    double ds = S0max/nbS;
    std::vector<double> S0_list(nbS + 1), prix_option(nbS +1);
    for (size_t k = 0; k <= nbS; ++k){
        S0_list[k] = k*ds;
        prix_option[k] =  prixAutocallableUnivariate(S0_list[k], Sref, B, r, sigma, T, nbT, Q_list, rbmt, Nmc); 
    }
    export_data(S0_list, {prix_option}, "graphAutocallableUnivariate");
}


double AutocallableUnivariate::deltaAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    std::vector<double> t(nbT), Sp(nbT), Sm(nbT);
    double phi_plus, phi_minus;
    double gauss;
    for (size_t k = 0; k < nbT; ++k)
        t[k] = (k+1);
    double estimateur = 0;
    for (size_t j = 0; j < Nmc; ++j){
        Sp.clear();
        Sm.clear();
        for (size_t i = 0; i < nbT; ++i){
            gauss = g();
            Sp[i] = (S0 + h) * exp( (r - sigma*sigma/2) * t[i] + sigma * sqrt(t[i]) * gauss );
            Sm[i]= (S0 - h) * exp( (r - sigma*sigma/2) * t[i] + sigma * sqrt(t[i]) * gauss );
        }
        phi_plus = verif(Sp, Sref, B, T, nbT, r, t, Q_list, rbmt);
        phi_minus = verif(Sm, Sref, B, T, nbT, r, t, Q_list, rbmt);
        estimateur += ( ( phi_plus - phi_minus ) / ( 2*h) );
    }
    estimateur /= Nmc;
    
    return estimateur;
}


void AutocallableUnivariate::graphDeltaAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    double ds = S0max/nbS;
    std::vector<double> S0_list(nbS + 1), deltaAutocallable_list(nbS +1);
    for (size_t k = 0; k <= nbS; ++k){
        S0_list[k] = k*ds;
        deltaAutocallable_list[k] =  deltaAutocallableUnivariate(S0_list[k], Sref, B, r, sigma, h, T, nbT, Q_list, rbmt, Nmc); 
    }
    export_data(S0_list, {deltaAutocallable_list}, "graphDeltaAutocallableUnivariate");
}


double AutocallableUnivariate::gammaAutocallableUnivariate(const double & S0, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    std::vector<double> t(nbT), Sp(nbT), Sm(nbT), S(nbT);
    double phi_plus, phi_minus, phi;
    double gauss;
    for (size_t k = 0; k < nbT; ++k)
        t[k] = (k+1);
    double estimateur = 0;
    for (size_t j = 0; j < Nmc; ++j){
        Sp.clear();
        Sm.clear();
        for (size_t i = 0; i < nbT; ++i){
            gauss = g();
            Sp[i] = (S0 + h) * exp( (r - sigma*sigma/2) * t[i] + sigma * sqrt(t[i]) * gauss );
            S[i] = S0 * exp( (r - sigma*sigma/2) * t[i] + sigma * sqrt(t[i]) * gauss );
            Sm[i]= (S0 - h) * exp( (r - sigma*sigma/2) * t[i] + sigma * sqrt(t[i]) * gauss );
        }
        phi_plus = verif(Sp, Sref, B, T, nbT, r, t, Q_list, rbmt);
        phi = verif(S, Sref, B, T, nbT, r, t, Q_list, rbmt);
        phi_minus = verif(Sm, Sref, B, T, nbT, r, t, Q_list, rbmt);
        estimateur += ( ( phi_plus - 2 * phi + phi_minus ) / ( h*h) );
    }
    estimateur /= Nmc;
    return estimateur;
}


void AutocallableUnivariate::graphGammaAutocallableUnivariate(const double & S0max, size_t const & nbS, const double & Sref, const double & B, const double & r, const double & sigma, const double & h, const int & T, size_t const & nbT, std::vector<double> const & Q_list, double (*rbmt)(double), const int & Nmc) const{
    double ds = S0max/nbS;
    std::vector<double> S0_list(nbS + 1), gammaAutocallable_list(nbS +1);
    for (size_t k = 0; k <= nbS; ++k){
        S0_list[k] = k*ds;
        gammaAutocallable_list[k] =  gammaAutocallableUnivariate(S0_list[k], Sref, B, r, sigma, h, T, nbT, Q_list, rbmt, Nmc); 
    }
    export_data(S0_list, {gammaAutocallable_list}, "graphGammaAutocallableUnivariate");
}


int main(){

    double S0(350), S0max(550), K(200), B(1), T(5), r(0.04), sigma(0.3), h(0.1), Sref(400);
    int Nmc(1000000), nbS(100), nbT(5);

    std::vector<double> Q_list = {110, 120, 130, 140, 150};
    
    AutocallableUnivariate univar;

    univar.graphAutocallableUnivariate(S0max, nbS, Sref, B, r, sigma, T, nbT, Q_list, rbmt1, Nmc);
    
    univar.graphDeltaAutocallableUnivariate(S0max, nbS, Sref, B, r, sigma, h, T, nbT, Q_list, rbmt1, Nmc);
    
    univar.graphGammaAutocallableUnivariate(S0max, nbS, Sref, B, r, sigma, h, T, nbT, Q_list, rbmt1, Nmc);

    return 0;
}