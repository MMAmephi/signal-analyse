#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <numbers>

using namespace std;

int sine_signal(int N_samples, double T, double phi, double F=10.0, int quantization_levels_num=256, double quantization_min = -1., double quantization_max = 1.){
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n";
    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((sin(2*3.1415*i/F/T + phi)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
    return 0;
}

int triangular_signal(int N_samples, double T, double phi, double t, double F=10.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n";
    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val=0;
        double* psv=&signal_val;
        if((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T < t){
            *psv = trunc((-abs((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;
        }
        else{
            *psv = 0;
        }
        csv_file << (i / F) << "," << *psv << "\n";
    }
    csv_file.close();
    return 0;
}

int rectangular_signal(int N_samples, double T, double phi, double t, double A, double F=10.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n";
    for (int i = 0; i < N_samples; ++i)
    {
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
        csv_file << (i / F) << "," << *psv << "\n";
    }
    csv_file.close();
    return 0;
}

int random_signal(int N_samples, double t, double A, double F=10.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n";
    double t1_start = (rand()%(N_samples/2-int(t*F)))/F;
    double t1_end = t1_start + t;
    double t2_start = (rand()%(N_samples-int(t*F))+trunc(t1_end*F))/F;
    double t2_end = t2_start + t;
    cout<< t1_start <<","<< t1_end <<","<< t2_start <<","<< t2_end <<"\n";
    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val=0;
        double* psv=&signal_val;
        if((i/F < t1_end && i/F >= t1_start) || (i/F < t2_end && i/F >= t2_start)){
            *psv = trunc((A)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;         
        }
        else{
            *psv = 0;
        }
        csv_file << (i / F) << "," << *psv << "\n";
    }
    csv_file.close();
    return 0;
}

int main()
{
    double F = 10.; // частота дискретизации
    int quantization_levels_num = 256; // количество уровней равномерного квантования
    double quantization_min = 0., quantization_max = 1.;
    int N_samples = 50;
    //int digital_signal[6] = {}; // коды квантования цифрового сигнала
    double phi = 0.; //Начальная фаза
    double T = 2.; //Период
    double t=0.5; //Длительность импульса
    double A=3.; //Амплитуда импульса
    //sine_signal(50, 2.0, 0.5);
    //triangular_signal(50, 2.0, 0.5, 1.0);
    rectangular_signal(128, 4.0, 3.14/2, 0.5, 3.0);
    //random_signal(100, 1.0, 3.0);
    return 0;
}