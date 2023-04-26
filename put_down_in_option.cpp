#include "utilities.hpp"


class DownAndIn{
    public:
    mutable NormalRandomGenerator g;
    double payOffPutEu(const double & S, const double & K) const;
    double payOffPutDownAndInEuS0fixe(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma) const;
    double prixPutDownAndInEuS0fixe(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc) const;  
    double payOffPutDownAndInEuSt(const double & B, const double & t, const double & St, const double & K, const double & T, const int & N, const double & r, const double & sigma) const; 
    double prixPutDownAndInEuSt(const double & B, const double & t, const double & St, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc) const;        
    void graphPayOffPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc ) const;
    void surfacePayOffPutDownAndIn(const double & B, const double & S0max, const int & nbS, const int & nbT, const double & K, const double & T, const int & N, const double & r,const double & sigma, const int & Nmc ) const;
    double deltaPut(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
    double deltaPutDownAndIn(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
    void graphDeltaPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
    double gammaPutDownAndIn(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
    double gammaPut(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
    void graphGammaPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const;
};


double DownAndIn::payOffPutEu(const double & S, const double & K) const{
    return std::max(K - S, 0.);
}


double DownAndIn::payOffPutDownAndInEuS0fixe(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma) const{
    if (S0 <= B)
            return 0;
    double St(S0), dt(T/N);
    bool flag = false;
    for (int i = 0; i < N; i++){
            St *= exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g());
            if (St <= B)
                flag = true;
    }
    if (flag)
        return payOffPutEu(St, K);
    return 0;
}


double DownAndIn::prixPutDownAndInEuS0fixe(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc) const{
    double gain(0), prix;
    gain = 0;
    for (int i = 0; i < Nmc; i++)
            gain += payOffPutDownAndInEuS0fixe(B, S0, K, T, N, r, sigma);
    prix = exp( -r * T ) * gain / Nmc;
    return prix;
}    


double DownAndIn::payOffPutDownAndInEuSt(const double & B, const double & t, const double & St, const double & K, const double & T, const int & N, const double & r, const double & sigma) const{
    if (St <= B)
            return 0;
    double dt;
    dt = ( T - t) / N;
    bool flag = false;
    double S(St);
    for (int i = 0; i < N; i++){
            S *= exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g());
            if (S <= B)
                flag = true;
    }
    if (flag)
        return payOffPutEu(S, K);
    return 0;
}


double DownAndIn::prixPutDownAndInEuSt(const double & B, const double & t, const double & St, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc) const{
    double gain, prix;
    gain = 0;
    for (int i = 0; i < Nmc; i++)
            gain += payOffPutDownAndInEuSt(B, t, St, K, T, N, r, sigma);
    prix = exp( -r * ( T - t) ) * gain / Nmc;
    return prix;
}


void DownAndIn::graphPayOffPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc) const{
    std::vector<double> S0_list(nbS + 1), pricePDI(nbS + 1)/*, prix(nbS + 1)*/;
    double ds = S0max/nbS;
    for (int k = 0; k <= nbS; k++){
            S0_list[k] = k*ds;
            pricePDI[k] = prixPutDownAndInEuS0fixe(B, S0_list[k], K, T, N, r, sigma, Nmc);
    }
    export_data(S0_list, {pricePDI}, "graphPayOffPutDownAndIn");
}


void DownAndIn::surfacePayOffPutDownAndIn(const double & B, const double & S0max, const int & nbS, const int & nbT, const double & K, const double & T, const int & N, const double & r, const double & sigma, const int & Nmc ) const{
    
    std::vector<double> t_list(nbT + 1), S0_list(nbS + 1);
    std::vector<std::vector<double>> prix(nbT + 1, std::vector<double>(nbS + 1));
    double dt, ds;
    dt = T/nbT;
    ds = S0max/nbS;
    for (int k = 0; k <= nbT; k++)
            t_list[k] = k*dt;        
    for (int k = 0; k <= nbS; k++)
            S0_list[k] = k*ds;
    for (int n = 0; n <= nbT; n++){
            for (int i = 0; i <= nbS; i++){
            prix[n][i] = prixPutDownAndInEuSt(B, t_list[n], S0_list[i], K, T, N, r, sigma, Nmc);
            }
    }
    export_data_surface(t_list, S0_list, prix, "surfacePayOffPutDownAndIn" );
    /*
    plt.plotSurface("$temps t$", "$prix actif S_{t}$", "$prix option V(0,S_{0})$", "Surface de put européen", t_list, nbT, S0_list, nbS, prix);
    const char* cmd1 = "python3 test.py surf";
    const int & result = system(cmd1);
    cout << "Result : " << result << std::endl;*/
}


double DownAndIn::deltaPutDownAndIn(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const{
    double Splus, Smoins, fact, dt, phi_plus, phi_moins, phi_res, delta_tmp,  delta;
    bool bp, bm;
    delta_tmp = 0;
    dt = T / N;
    for (int k = 0; k < Nmc; k++){
            Splus = S0 + h;
            Smoins = S0 - h;
            bp = bm = false;
            phi_plus = phi_moins = 0;
            if (Smoins > B){
            for (int j = 0; j < N; j++){
                    fact = exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g() );
                    Splus *= fact;
                    Smoins *= fact;
                    if (Splus <= B)
                    bp = true;
                    if (Smoins <= B)
                    bm = true;
            }
            if (bp)
                    phi_plus = payOffPutEu(Splus, K);
            if (bm)
                    phi_moins = payOffPutEu(Smoins, K);
            }
            else {
            if (Splus > B){
                    for (int j = 0; j < N; j++){
                    fact = exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g() );
                    Splus *= fact;
                    if (Splus <= B)
                            bp = true;
                    }
                    if (bp)
                    phi_plus = payOffPutEu(Splus, K);
            }
            }
            phi_res = ( phi_plus - phi_moins) / ( 2 * h );
            delta_tmp += phi_res;
    }
    delta = exp( -r * T ) * delta_tmp / Nmc;
    return delta;
}


void DownAndIn::graphDeltaPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS + 1), delta_list(nbS + 1), delta(nbS + 1);
    double ds = S0max/nbS;
    for (int k = 0; k <= nbS; k++){
            S0_list[k] = k*ds;
            delta_list[k] = deltaPutDownAndIn(B, S0_list[k], K, T, N, r, sigma, h, Nmc);
    }
    export_data(S0_list, {delta_list},"graphDeltaPutDownAndIn");
}


double DownAndIn::gammaPutDownAndIn(const double & B, const double & S0, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const{
    double Splus, S, Smoins, fact, dt, phi_plus, phi, phi_moins, phi_res, gamma_tmp,  gamma;
    bool bp, b, bm;
    gamma_tmp = 0;
    dt = T / N;
    for (int k = 0; k < Nmc; k++){
            Splus = S0 + h;
            S = S0;
            Smoins = S0 - h;
            bp = b = bm = false;
            phi_plus = phi_moins = 0;
            if (Smoins > B){
            for (int j = 0; j < N; j++){
                    fact = exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g() );
                    Splus *= fact;
                    S *= fact;
                    Smoins *= fact;
                    if (Splus <= B)
                    bp = true;
                    if (S <= B)
                    S = true;
                    if (Smoins <= B)
                    bm = true;
            }
            if (bp)
                    phi_plus = payOffPutEu(Splus, K);
            if (b)
                    phi = payOffPutEu(S, K);
            if (bm)
                    phi_moins = payOffPutEu(Smoins, K);
            }
            else {
            if (S > B){
                    for (int j = 0; j < N; j++){
                    fact = exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g() );
                    Splus *= fact;
                    S *= fact;

                    if (Splus <= B)
                            bp = true;
                    if (S <= B)
                            S = true;
                    }
                    if (bp)
                    phi_plus = payOffPutEu(Splus, K);
                    if (b)
                    phi = payOffPutEu(S, K);
            }
            else {
                    if (Splus > B){
                    for (int j = 0; j < N; j++){
                            fact = exp( ( r - sigma*sigma / 2 ) * dt + sigma * sqrt(dt) * g() );
                            Splus *= fact;
                            if (Splus <= B)
                            bp = true;
                    }
                    if (bp)
                            phi_plus = payOffPutEu(Splus, K);
                    }
            }
        }
        phi_res = ( phi_plus - 2 * phi + phi_moins) / ( h * h );
        gamma_tmp += phi_res;
    }
    gamma = exp( -r * T ) * gamma_tmp / Nmc;
    return gamma;
}


void DownAndIn::graphGammaPutDownAndIn(const double & B, const double & S0max, const int & nbS, const double & K, const double & T, const int & N, const double & r, const double & sigma, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS + 1), gamma_list(nbS + 1);
    double ds = S0max/nbS;
    for (int k = 0; k <= nbS; k++){
            S0_list[k] = k*ds;
            gamma_list[k] = gammaPutDownAndIn(B, S0_list[k], K, T, N, r, sigma, h, Nmc);
    }
    export_data(S0_list, {gamma_list},"graphGammaPutDownAndIn");
}


int main(){

    double S0(100), S0max(550), K(200), B(80), T(5), r(0.04), sigma(0.3), h(0.1), Sref(400), tt(3), Stt(100);
    int Nmc(1000000), nbS(100), nbT(5), N(100);
    std::vector<double> Q_list = {110, 120, 130, 140, 150};
    std::vector<double> t = {1, 2, 3, 4, 5};

    DownAndIn dAi;

    dAi.graphPayOffPutDownAndIn(B, S0max, nbS, K, T, N, r, sigma, Nmc);
    
    dAi.surfacePayOffPutDownAndIn(B, S0max, nbS, nbT, K, T, N, r, sigma, Nmc);    
    
    dAi.graphDeltaPutDownAndIn(B, S0max, nbS, K, T, N, r, sigma, h, Nmc);
    
    dAi.graphGammaPutDownAndIn(B, S0max, nbS, K, T, N, r, sigma, h, Nmc);
    
    return 0;
}