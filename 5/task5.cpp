#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <complex.h>
#include <vector>
#include <random>

using namespace std;

void vector_csv(vector<double> data, int N, double F){
    int N_samples = (int) N*F;
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "t, signal\n";
    for(int i=0; i <= N_samples; i++){
        csv_file << (i / F) << "," << data[i] << "\n";    
    }    
    csv_file.close();
}

vector<double> sine_signal(int N_samples, double T, double phi, double F=16.0, int quantization_levels_num=256, double quantization_min = -1., double quantization_max = 1.){
    vector<double> signal;
    for (int i = 0; i <= N_samples; ++i){
        double signal_val = trunc((sin(2*3.1415*i/F/T + phi)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        signal.push_back(signal_val);
    }
    return signal;    
}

vector<double> triangular_signal(int N_samples, double T, double phi, double t, double F=16.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    vector<double> signal;
    for (int i = 0; i <= N_samples; ++i){
        double signal_val=0;
        double* psv=&signal_val;
        if((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T < t){
            *psv = trunc((-abs((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;
        }
        else{
            *psv = 0;
        }
        signal.push_back(*psv);
    }
    return signal;
}

vector<double> rectangular_signal(int N_samples, double T, double phi, double t, double A, double F=16.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    vector<double> signal;
    for (int i = 0; i <= N_samples; ++i){
        double signal_val=0;
        double* psv=&signal_val;
        if(i==0){
            *psv = 0;
        }
        else if((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T < t){
            *psv = trunc((A)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;           
        }
        else{
            *psv = 0;
        }
        signal.push_back(*psv);
    }
    return signal;
}

vector<double> random_signal(int N_samples, double t, double A, double F=16.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    vector<double> signal;
    double t1_start = (rand()%(N_samples/2-int(t*F)))/F;
    double t1_end = t1_start + t;
    double t2_start = (rand()%(N_samples-int(t*F))+trunc(t1_end*F))/F;
    double t2_end = t2_start + t;
    //cout<< t1_start <<","<< t1_end <<","<< t2_start <<","<< t2_end <<"\n";
    for (int i = 0; i <= N_samples; ++i)
    {
        double signal_val=0;
        double* psv=&signal_val;
        if((i/F < t1_end && i/F >= t1_start) || (i/F < t2_end && i/F >= t2_start)){
            *psv = trunc((A)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;         
        }
        else{
            *psv = 0;
        }
        signal.push_back(*psv);
    }
    return signal;
}

vector<double> step_interpolation(vector<double> signal, int N_samples, double a, double b, int N_samples_grid){
    vector<double> signal_grid;
    double h = (b - a) / N_samples; 
    double h_grid = (b - a) / N_samples_grid;
    double temp; 
    int pos;
    for(int i = 0; i <= N_samples_grid; i++){
        pos = trunc((i*h_grid)/h);
        temp = signal[pos];
        signal_grid.push_back(temp);
    }
    return signal_grid;
}

vector<double> linear_interpolation(vector<double> signal, int N_samples, double a, double b, int N_samples_grid){
    vector<double> signal_grid;
    double h = (b - a) / N_samples; 
    double h_grid = (b - a) / N_samples_grid;
    double temp;
    int pos; 
    for(int i = 0; i < N_samples_grid; i++){
        pos = trunc((i*h_grid)/h);
        temp = (signal[pos] + (signal[pos+1]-signal[pos])/h*(i*h_grid - pos*h));
        signal_grid.push_back(temp);
    }
    return signal_grid;
}

int main() {
    double F = 4.;
    double N = 8;
    int N_samples = (int) N * F;
    double signal_value;
    double phi = 3.14 / 2.;
    double period = 2.;
    double t=0.5; //Длительность импульса
    double A=3.;

    double a = -1.;
    double b = 1.;

    double F_grid = 16.;
    int N_samples_grid =(int) N * F_grid;

    vector<double> temp = sine_signal(N_samples, period, phi, F=F);
    //vector<double> temp = rectangular_signal(N_samples, period, phi, t, A);
    //vector<double> temp = triangular_signal(N_samples, period, phi, t);
    //vector<double> temp = random_signal(N_samples, t, A);

    vector<double> signal_grid;

    //signal_grid = step_interpolation(temp, N_samples, a, b, N_samples_grid);
    signal_grid = linear_interpolation(temp, N_samples, a, b, N_samples_grid);

    //vector_csv(temp, N, F);
    vector_csv(signal_grid, N, F_grid);
    return 0;
}