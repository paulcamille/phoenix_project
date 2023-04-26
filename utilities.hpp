#include <iostream> 
#include <vector>
#include <random>
#include <fstream>  
#include <algorithm>
#include <iterator>  


class NormalRandomGenerator {
    public:
        double operator()();

    private:
        static std::mt19937 generator;
        static std::normal_distribution<double> distribution;
};


std::mt19937 NormalRandomGenerator::generator(std::random_device{}());
std::normal_distribution<double> NormalRandomGenerator::distribution(0.0, 1.0);


double NormalRandomGenerator::operator()() {
    return distribution(generator);
}


void export_data(const std::vector<double>  & t1, const std::vector<double>  & t2, std::string nameGraph){
    std::ofstream f_data;
    f_data.open("data/"+nameGraph+".csv");
    if (f_data){
        for (int n = 0; n < t1.size(); n++)
            f_data << t1[n] << "\t" << t2[n] << std::endl;
        f_data.close();
    }
    else
        std::cout << "erreur" << std::endl;  
}


void export_data_surface(std::vector<double> const & t1, std::vector<double> const &  t2, std::vector<std::vector<double>> const & t3, std::string nameGraph){

    std::ofstream f_data;
    f_data.open("data/"+nameGraph+".csv");
    size_t l1 = t1.size();
    size_t l2= t2.size();
    if (f_data){
        size_t n, m;
        for (m = 0; m < l2; m++)
            f_data << "\t" << t2[m];
        f_data << std::endl;
        for (n = 0; n < l1; ++n){
            f_data << t1[n];
            for (m = 0; m < l2; ++m)
                f_data << "\t" << t3[n][m];
            f_data << std::endl;
        }
        f_data.close();
    }
    else
        std::cout << "erreur" << std::endl;  
}