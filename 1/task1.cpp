#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <numbers>

using namespace std;

int main()
{
    double F = 10.; // частота дискретизации
    int quantization_levels_num = 256; // количество уровней равномерного квантования
    double quantization_min = 0., quantization_max = 1.;
    int N_samples = 50;
    //int digital_signal[6] = {}; // коды квантования цифрового сигнала
    const char csv_file_name[64] = "data.csv"; // имя файла для вывода

    double phi = 0.; //Начальная фаза
    double T = 2.; //Период
    double t=0.5; //Длительность импульса
    double A=3.; //Амплитуда импульса

    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n"; // записать заголовки колонок
    /*for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((sin(2*M_PI*i/F/T + phi)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }*/
    /*for (int i = 0; i < N_samples; ++i)
    {
        double signal_val=0;
        double* psv=&signal_val;
        if((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T < t){
            *psv = trunc((-abs((i/F/T+phi/T)*T-trunc(i/F/T+phi/T)*T-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;
        }
        else{
            *psv = 0;
        }
        //cout <<  trunc((-abs(i/F-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) << ", " << signal_val << "\n";
        csv_file << (i / F) << "," << *psv << "\n";
    }*/
    /*for (int i = 0; i < N_samples; ++i)
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
        //cout <<  trunc((-abs(i/F-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) << ", " << signal_val << "\n";
        csv_file << (i / F) << "," << *psv << "\n";
    }*/
    double t1_start = (rand()%(N_samples/2-int(t*F)))/F;
    double t1_end = t1_start + t;
    double r2_start = (rand()%N_samples+)
    cout<< t1_start <<","<< t1_end <<","<< t2_start <<","<< t2_end <<",";
    csv_file.close();
    return 0;
}