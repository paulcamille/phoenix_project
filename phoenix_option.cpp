#include "utilities.hpp"


class Phoenix{
    public:
        mutable NormalRandomGenerator g;
        Phoenix();

        double payOffPhoenix(const std::vector<double>  & S0walk, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const;
        double getPayOffPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const;
        double pricePhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const;
        void graphPricePhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const;    
        void graphPricePhoenixProportionalBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const;
        void graphPricePhoenixCYnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph,  const double & CYmax, const int & nbCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixCphnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & Cphmax, const int & nbCph,const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixBphnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & Bphmax,  const int & nbBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixBYnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & BYmax, const int & nbBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixBputnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & Bputmax,  const int & nbBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixSigmanotFixS0fix(const double & S0, const double & r, const double & sigmamax, const int & nbSigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const;
        void graphPricePhoenixRnotFixS0fix(const double & S0, const double & rmax,  const int & nbR, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const;
        void graphPricePhoenixTnotFixS0fix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & Tmax, const int & Nmc) const;
        double deltaPhoenix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        void graphDeltaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        double gammaPhoenix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        void graphGammaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        void partition(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc, const int N, const double & alpha) const;
        
        double payOffPhoenixNotActualised(const std::vector<double>  & S0walk, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const;
        double pricePhoenixNotActualised(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const;
        double getPayOffNotActualised(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const;

        double vegaPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        void graphVegaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;
        void graphVegaPhoenixSigmanotFixS0fix(const double & S0, const double & r, const double & sigmamax, const int & nbSigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc)const;
        void graphVegaPhoenixProportionalBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const;




};

Phoenix::Phoenix() {}


double Phoenix::payOffPhoenix(const std::vector<double>  & S0walk, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const{ 
    bool realFlag = true;
    double payOff = 0;
    size_t i = 1;
    while (i < nbT && realFlag){
        if (S0walk[i] > Bph){
            realFlag = false;
            payOff += (No + Cph) * exp(-r*(t[i]-t[0]));
        }
        else{
            if(S0walk[i] > BY)
                payOff += CY * exp(-r*(t[i]-t[0]));
        }
        ++i;
    }
    if (realFlag && (i == nbT)){
        if (S0walk[i] > Bph){
            realFlag = false;
            payOff += (No + Cph) * exp(-r*(T-t[0]));
        }
        else{
            if(S0walk[i] > BY)
                payOff += (No + CY) * exp(-r*(T-t[0]));
            else{
                if(S0walk[i] > Bput)
                    payOff += No * exp(-r*(T-t[0]));
                else
                    payOff += std::max((K - S0walk[i])/S0walk[0], 0.0) * exp(-r*(T-t[0]));
            }
        }
    }
    return payOff;
}


double Phoenix::payOffPhoenixNotActualised(const std::vector<double> & S0walk, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double> & t, const int & nbT) const{ 
    bool realFlag = true;
    double payOff = 0;
    size_t i = 1;
    while (i < nbT && realFlag){
        if (S0walk[i] > Bph){
            realFlag = false;
            payOff += (No + Cph) * exp(-r*(t[i]- T));
        }
        else{
            if(S0walk[i] > BY)
                payOff += CY * exp(-r*(t[i]- T));
        }
        ++i;
    }
    if (realFlag && (i == nbT)){
        if (S0walk[i] > Bph){
            realFlag = false;
            payOff += (No + Cph);
        }
        else{
            if(S0walk[i] > BY)
                payOff += (No + CY) ;
            else{
                if(S0walk[i] > Bput)
                    payOff += No ;
                else
                    payOff += std::max((K - S0walk[i])/S0walk[0], 0.0) ;
            }
        }
    }
    return payOff;
}


double Phoenix::getPayOffPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const{
    std::vector<double> S0walk(nbT+1);
    S0walk[0] = S0;
    for (size_t k = 1; k <= nbT; ++k)
        S0walk[k] = S0 * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * g());
    return payOffPhoenix(S0walk, r, sigma, K, BY, Bput, Bph, CY, Cph, No, T, t, nbT);
}


double Phoenix::getPayOffNotActualised(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT) const{
    std::vector<double> S0walk(nbT+1);
    S0walk[0] = S0;
    for (size_t k = 1; k <= nbT; ++k)
        S0walk[k] = S0 * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * g());
    return payOffPhoenixNotActualised(S0walk, r, sigma, K, BY, Bput, Bph, CY, Cph, No, T, t, nbT);
}


double Phoenix::pricePhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const{
    double price = 0;
    std::vector<double> S0walk(nbT+1); 
    for (size_t i = 0; i < Nmc; ++i){
        price += getPayOffPhoenix(S0, r, sigma, K, BY, Bput, Bph, CY, Cph, No, T, t, nbT);
    }
    return price/Nmc;
}


double Phoenix::pricePhoenixNotActualised(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const{ 
    double price = 0;
    std::vector<double> S0walk(nbT+1); 
    for (size_t i = 0; i < Nmc; ++i){
        price += getPayOffNotActualised(S0, r, sigma, K, BY, Bput, Bph, CY, Cph, No, T, t, nbT);
    }
    return price/Nmc;
}


void Phoenix::graphPricePhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const{
    double ds = S0max/nbS;
    std::vector<double> S0_list(nbS), price(nbS);
    for (size_t k = 0; k < nbS; ++k){
        S0_list[k] = (k+1)*ds;
        price[k] =  pricePhoenix(S0_list[k], r, sigma, S0_list[k], BY, Bput, Bph, CY, Cph, No, T, t, nbT, Nmc);
    }
    export_data(S0_list, price, "graphPricePhoenixFixBarrier");
}


void Phoenix::graphPricePhoenixProportionalBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const{
    double ds = S0max/nbS;
    std::vector<double> S0_list(nbS), price(nbS);
    for (size_t k = 0; k < nbS; ++k){
        S0_list[k] = (k+1)*ds;
        price[k] =  pricePhoenix(S0_list[k], r, sigma, S0_list[k]*coefK, S0_list[k]*coefBY, S0_list[k]*coefBput, S0_list[k]*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(S0_list, price, "graphPricePhoenixProportionalBarrier");
}


void Phoenix::graphPricePhoenixCYnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph,  const double & CYmax, const int & nbCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dcy = CYmax/nbCY;
    std::vector<double> CY_list(nbCY+1), price(nbCY+1);
    for (size_t k = 0; k <= nbCY; ++k){
        CY_list[k] = k*dcy;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*coefBY, S0*coefBput, S0*coefBph, No*CY_list[k], No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(CY_list, price, "graphPricePhoenixCYnotFix");
}


void Phoenix::graphPricePhoenixCphnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & Cphmax, const int & nbCph ,const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dcph = Cphmax/nbCph;
    std::vector<double> Cph_list(nbCph+1), price(nbCph+1);
    for (size_t k = 0; k <= nbCph; ++k){
        Cph_list[k] = (k+1)*dcph;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*coefBY, S0*coefBput, S0*coefBph, No*coefCY, No*Cph_list[k], No, T, t, nbT, Nmc);
    }
    export_data(Cph_list, price, "graphPricePhoenixCphnotFix");
}


void Phoenix::graphPricePhoenixBphnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & Bphmax,  const int & nbBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dcph = Bphmax/nbBph;
    std::vector<double> Bph_list(nbBph+1), price(nbBph+1);
    for (size_t k = 0; k <= nbBph; ++k){
        Bph_list[k] = (k+1)*dcph;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*coefBY, S0*coefBput, S0*Bph_list[k], No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(Bph_list, price, "graphPricePhoenixBphnotFix");
}


void Phoenix::graphPricePhoenixBYnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & BYmax, const int & nbBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dcph = BYmax/nbBY;
    std::vector<double> BY_list(nbBY+1), price(nbBY+1);
    for (size_t k = 0; k <= nbBY; ++k){
        BY_list[k] = (k+1)*dcph;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*BY_list[k], S0*coefBput, S0*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(BY_list, price, "graphPricePhoenixBYnotFix");
}


void Phoenix::graphPricePhoenixBputnotFix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & Bputmax,  const int & nbBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dcph = Bputmax/nbBput;
    std::vector<double> Bput_list(nbBput+1), price(nbBput+1);
    for (size_t k = 0; k <= nbBput; ++k){
        Bput_list[k] = (k+1)*dcph;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*coefBY, S0*Bput_list[k], S0*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(Bput_list, price, "graphPricePhoenixBputnotFix");
}


void Phoenix::graphPricePhoenixSigmanotFixS0fix(const double & S0, const double & r, const double & sigmamax, const int & nbSigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc)const{
    double dsigma = sigmamax/nbSigma;
    std::vector<double> Sigma_list(nbSigma), price(nbSigma);
    for (size_t k = 0; k < nbSigma; ++k){
        Sigma_list[k] = (k+1)*dsigma;
        price[k] =  pricePhoenix(S0, r, Sigma_list[k], S0*coefK, S0*coefBY, S0*coefBput, S0*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(Sigma_list, price, "graphPricePhoenixSigmanotFixS0fix");
}


void Phoenix::graphPricePhoenixRnotFixS0fix(const double & S0, const double & rmax,  const int & nbR, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const int & Nmc) const{
    double dr = rmax/nbR;
    std::vector<double> R_list(nbR), price(nbR);
    for (size_t k = 0; k < nbR; ++k){
        R_list[k] = (k+1)*dr;
        price[k] =  pricePhoenix(S0,  R_list[k], sigma, S0*coefK, S0*coefBY, S0*coefBput, S0*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, Nmc);
    }
    export_data(R_list, price, "graphPricePhoenixRnotFixS0fix");
}


void Phoenix::graphPricePhoenixTnotFixS0fix(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & Tmax, const int & Nmc) const{
    double dt = Tmax/Tmax;
    std::vector<double> T_list(Tmax), price(Tmax);
    auto tab = [&] (const int & ed){
        std::vector<double> timelist(ed+1);
        for (size_t k = 0; k <= ed; ++k)
            timelist[k] = k;
        return timelist;
    };
    for (size_t k = 0; k < Tmax; ++k){
        T_list[k] = (k+1)*dt;
        price[k] =  pricePhoenix(S0, r, sigma, S0*coefK, S0*coefBY, S0*coefBput, S0*coefBph, No*coefCY, No*coefCph, No, T_list[k], tab(T_list[k]), T_list[k], Nmc);
    }
    export_data(T_list, price, "graphPricePhoenixTnotFixS0fix");
}


double Phoenix::deltaPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0walkplus(nbT+1), S0walkminus(nbT+1);
    double delta(0), s_p(S0+h), s_m(S0-h);
    double gauss;
    for (int k = 0; k < Nmc; k++){
        S0walkplus[0] = s_p;
        S0walkminus[0] = s_m;
        for (size_t k = 1; k <= nbT; ++k){
            gauss = g();
            S0walkplus[k] = s_p * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * gauss);
            S0walkminus[k] = s_m * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * gauss);
        }
        delta += (payOffPhoenix(S0walkplus, r, sigma, s_p, BY, Bput, Bph, CY, Cph, No, T, t, nbT) - payOffPhoenix(S0walkminus, r, sigma, s_m, BY, Bput, Bph, CY, Cph, No, T, t, nbT)) / (2*h);
    }
    return delta/Nmc;
}


void Phoenix::graphDeltaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS), delta_list(nbS);
    double ds = S0max/nbS;
    for (int k = 0; k < nbS; k++){
        S0_list[k] = (k+1)*ds;
        delta_list[k] = deltaPhoenix(S0_list[k], r, sigma, K, BY, Bput, Bph, CY, Cph, No, T, t, nbT, h, Nmc);
    }
    export_data(S0_list, delta_list,"graphDeltaPhoenixFixBarrier");
}


double Phoenix::gammaPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double> & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0walkplus(nbT+1), S0walkminus(nbT+1), S0walk(nbT+1);
    double gamma(0), s_p(S0+h), s_m(S0-h);
    double gauss;
    for (int k = 0; k < Nmc; k++){
        S0walkplus[0] = (S0+h);
        S0walkminus[0] = S0-h;
        S0walk[0] = S0;
        for (size_t k = 1; k <= nbT; ++k){
            gauss = g();
            S0walkplus[k] = s_p * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * gauss);
            S0walkminus[k] = s_m * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * gauss);
            S0walk[k] = S0 * exp( (r - sigma*sigma/2) * (t[k]-t[0])+ sigma * sqrt(t[k]-t[0]) * gauss);
        }
        gamma += (payOffPhoenix(S0walkplus, r, sigma, s_p, BY, Bput, Bph, CY, Cph, No, T, t, nbT) + payOffPhoenix(S0walkminus, r, sigma, s_m, BY, Bput, Bph, CY, Cph, No, T, t, nbT) - 2 * payOffPhoenix(S0walk, r, sigma, S0, BY, Bput, Bph, CY, Cph, No, T, t, nbT)) / (h*h);
    }
    return gamma/Nmc;
}


void Phoenix::graphGammaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS), gamma_list(nbS);
    double ds = S0max/nbS;
    for (int k = 0; k < nbS; k++){
        S0_list[k] = (k+1)*ds;
        gamma_list[k] = gammaPhoenix(S0_list[k], r, sigma, S0_list[k], BY, Bput, Bph, CY, Cph, No, T, t, nbT, h, Nmc);
    }
    export_data(S0_list, gamma_list,"graphGammaPhoenixFixBarrier");
}


double Phoenix::vegaPhoenix(const double & S0, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0walkplus(nbT+1), S0walkminus(nbT+1);
    double vega(0), sigma_p(sigma+h), sigma_m(sigma-h);
    double gauss;
    for (int k = 0; k < Nmc; k++){
        S0walkplus[0] = S0;
        S0walkminus[0] = S0;
        for (size_t k = 1; k <= nbT; ++k){
            gauss = g();
            S0walkplus[k] = S0 * exp( (r - sigma_p*sigma_p/2) * (t[k]-t[0])+ sigma_p * sqrt(t[k]-t[0]) * gauss);
            S0walkminus[k] = S0 * exp( (r - sigma_m*sigma_m/2) * (t[k]-t[0])+ sigma_m * sqrt(t[k]-t[0]) * gauss);
        }
        vega += (payOffPhoenix(S0walkplus, r, sigma_p, S0, BY, Bput, Bph, CY, Cph, No, T, t, nbT) - payOffPhoenix(S0walkminus, r, sigma_m, S0, BY, Bput, Bph, CY, Cph, No, T, t, nbT)) / (2*h);
    }
    return vega/Nmc;
}


void Phoenix::graphVegaPhoenixFixBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & K, const double & BY, const double & Bput, const double & Bph, const double & CY, const double & Cph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS), vega_list(nbS);
    double ds = S0max/nbS;
    for (int k = 0; k < nbS; k++){
        S0_list[k] = (k+1)*ds;
        vega_list[k] = vegaPhoenix(S0_list[k], r, sigma, S0_list[k], BY, Bput, Bph, CY, Cph, No, T, t, nbT, h, Nmc);
    }
    export_data(S0_list, vega_list,"graphVegaPhoenixFixBarrier");
}


void Phoenix::graphVegaPhoenixSigmanotFixS0fix(const double & S0, const double & r, const double & sigmamax, const int & nbSigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc)const{
    double dsigma = sigmamax/nbSigma;
    double k_s0(coefK*S0), by_s0(coefBY*S0), bput_s0(coefBput*S0), bph_s0(coefBph*S0), cy_No(coefCY*No), cph_No(coefCph*No);
    std::vector<double> Sigma_list(nbSigma), vega_list(nbSigma);
    for (size_t k = 0; k < nbSigma; ++k){
        Sigma_list[k] = (k+1)*dsigma;
        vega_list[k] =  vegaPhoenix(S0, r, Sigma_list[k], k_s0, by_s0, bput_s0, bph_s0, cy_No, cph_No, No, T, t, nbT, h, Nmc);
    }
    export_data(Sigma_list, vega_list, "graphVegaPhoenixSigmanotFixS0fix");
}


void Phoenix::graphVegaPhoenixProportionalBarrier(const double & S0max, const int & nbS, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput,const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double>  & t, const int & nbT, const double & h, const int & Nmc) const{
    std::vector<double> S0_list(nbS), vega_list(nbS);
    double ds = S0max/nbS;
    for (int k = 0; k < nbS; k++){
        S0_list[k] = (k+1)*ds;
        vega_list[k] = vegaPhoenix(S0_list[k], r, sigma, S0_list[k]*coefK, S0_list[k]*coefBY, S0_list[k]*coefBput, S0_list[k]*coefBph, No*coefCY, No*coefCph, No, T, t, nbT, h, Nmc);
    }
    export_data(S0_list, vega_list,"graphVegaPhoenixProportionalBarrier");
}


std::vector<double> densite(std::vector<double> & tab){
    const int & ind = tab.size();
    std::sort(tab.begin(), tab.end());
    std::vector<double> distribution(ind);
    double const step = 1.0/ind;
    for (int i = 1; i <= ind; ++i){
        distribution[i-1] = i*step;
    }
    export_data(tab, {{distribution}}, "distribution");
    return std::move(distribution);
}


void Phoenix::partition(const double & S0, const double & r, const double & sigma, const double & coefK, const double & coefBY, const double & coefBput, const double & coefBph, const double & coefCY, const double & coefCph, const double & No, const int & T, const std::vector<double> & t, const int & nbT, const int & Nmc, const int N, const double & alpha) const{
    std::vector<double> vect(N);
    double k_s0(coefK*S0), by_s0(coefBY*S0), bput_s0(coefBput*S0), bph_s0(coefBph*S0), cy_No(coefCY*No), cph_No(coefCph*No);
    double V0(pricePhoenixNotActualised(S0, r, sigma, k_s0, by_s0, bput_s0, bph_s0, cy_No, cph_No, No, T, t, nbT, Nmc));
    for(size_t i = 0; i < N; ++i)
        vect[i] = getPayOffNotActualised(S0,r,sigma, k_s0, by_s0, bput_s0, bph_s0, cy_No, cph_No, No, T,  t, nbT) - V0;
    std::vector<double> distribution(densite(vect));
    double VaR(vect[int(alpha * N)-1]);
    std::ofstream write_on("data/var.csv");
    if (write_on){
        write_on << VaR << "\t" << distribution[int(alpha * N)-1] << std::endl;
        write_on.close();
    }
    else 
        std::cout << "error" << std::endl;
}


int main(){

    std::vector<double> t = {0, 1, 2, 3, 4, 5};

    double S0(100), S0max(200) ,r(0.02) ,sigma(0.3) ,K(S0) ,BY(80) ,Bput(70) ,Bph(120) , No(1), CY(0.05*No) ,Cph(0.1*No), coefK(100.0/100) ,coefBY(0.8),coefBput(0.7), coefBph(1.2), coefCY(0.05), coefCph(0.1), CYmax(5), Cphmax(5), Bphmax(5), BYmax(5), Bputmax(5), sigmamax(7.5), rmax(5), h(0.1), alpha(1.0/100);

    int nbS(100), T(5), N(100000), nbT(5), Nmc(1000000), nbCY(nbS), nbCph(nbS), nbBph(nbS), nbBY(nbS), nbBput(nbS), nbSigma(nbS), nbR(nbS), Tmax(20);
  
    Phoenix feu;

    std::cout <<  "payoff =\t" << feu.getPayOffPhoenix(S0,r,sigma, K ,BY,Bput, Bph, CY, Cph, No, T,  t, nbT) << std::endl;
 
    std::cout << "price =\t" << feu.pricePhoenix(S0,r,sigma,K ,BY,Bput, Bph, CY, Cph, No, T,  t, nbT, Nmc) << std::endl;;

    feu.graphPricePhoenixFixBarrier(S0max, nbS, r, sigma, K,  BY, Bput, Bph, CY, Cph, No, T,  t, nbT, Nmc);
    
    feu.graphPricePhoenixProportionalBarrier(S0max, nbS, r, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixCYnotFix(S0, r, sigma, coefK,  coefBY, coefBput, coefBph, CYmax, nbCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixCphnotFix(S0, r, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, Cphmax, nbCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixBphnotFix(S0, r, sigma, coefK,  coefBY, coefBput, Bphmax, nbBph, coefCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixBYnotFix(S0, r, sigma, coefK,  BYmax, nbBY, coefBput, Bph, coefCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixBputnotFix(S0, r, sigma, coefK,  coefBY, Bputmax, nbBput, coefBph, coefCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixSigmanotFixS0fix(S0, r, sigmamax, nbSigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, Nmc);

    feu.graphPricePhoenixRnotFixS0fix(S0, rmax, nbR, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, Nmc);
    
    feu.graphPricePhoenixTnotFixS0fix(S0, r, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, Tmax, Nmc);

    std::cout << "delta =\t" << feu.deltaPhoenix(S0, r, sigma, K ,BY,Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc) <<  std::endl;
    
    feu.graphDeltaPhoenixFixBarrier(S0max, nbS, r, sigma, K,  BY, Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc);

    std::cout << "gamma =\t" << feu.gammaPhoenix(S0, r, sigma, K ,BY,Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc) <<  std::endl;

    feu.graphGammaPhoenixFixBarrier(S0max, nbS, r, sigma, K,  BY, Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc);

    std::cout << "vega =\t" << feu.vegaPhoenix(S0, r, sigma, K ,BY,Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc) <<  std::endl;

    feu.graphVegaPhoenixFixBarrier(S0max, nbS, r, sigma, K,  BY, Bput, Bph, CY, Cph, No, T,  t, nbT, h, Nmc);

    feu.partition(S0, r, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, Nmc, N, alpha);

    feu.graphVegaPhoenixSigmanotFixS0fix(S0, r, sigmamax, nbSigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, h,Nmc);

    feu.graphVegaPhoenixProportionalBarrier(S0max, nbS, r, sigma, coefK,  coefBY, coefBput, coefBph, coefCY, coefCph, No, T,  t, nbT, h, Nmc);

    return 0;
}