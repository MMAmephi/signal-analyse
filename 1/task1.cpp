#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <numbers>

int main()
{
    double F = 10.; // частота дискретизации
    int quantization_levels_num = 64; // количество уровней равномерного квантования
    double quantization_min = 0., quantization_max = 1.;
    int N_samples = 50;
    //int digital_signal[6] = {}; // коды квантования цифрового сигнала
    const char csv_file_name[64] = "data.csv"; // имя файла для вывода

    double phi = 1.; //Начальная фаза
    double T = 2.; //Период
    double t=2.; //Длительность импульса

    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time, signal\n"; // записать заголовки колонок
    /*for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((sin(2*M_PI*i/F/T + phi)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }*/
    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((-abs(i/F-t/2)+t/2)*(quantization_levels_num - 1)/(quantization_max - quantization_min)) * ((quantization_max - quantization_min)/(quantization_levels_num - 1))+quantization_min;
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
    return 0;
}