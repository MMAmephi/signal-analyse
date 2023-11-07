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
    csv_file << "signal\n";
    for(int i=0; i < N_samples; i++){
        csv_file << (i / F) << "," << data[i] << "\n";    
    }    
    csv_file.close();
}

void vector_csv_complex(vector<complex<double>> data, int N, double F){
    int N_samples = (int) N*F;
    const char csv_file_name[64] = "data.csv";
    ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "signal\n";
    for(int i=0; i < N_samples; i++){
        csv_file << -F/2 + double(i)/N << "," << abs(data[i]) << "\n";    
    }    
    csv_file.close();
}

void AWGN_signal(vector<double> &signal, double s, double a, double b) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> X(0, s);
    double noise;
	for (int i = 0; i < signal.size(); i++) {
        noise = X(gen);
        if(((signal[i]+noise) < b) && ((signal[i]+noise) > a)){
            signal[i] += noise;
        }
	}
}

void Impulse_signal(vector<double> &signal, double P, double a, double b) {
	double tmp;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distrib(0, 100);
    for (int i = 0; i < signal.size(); i++) {
        tmp = distrib(gen);
        if (tmp <= (100*P / 2)) {
		    signal[i] = a;
		}
		if ((tmp <= (100*P)) && (tmp > (100 * P / 2))) {
			signal[i] = b;
		}
	}
	return;
}

vector<int> hist(vector<double> signal, int bins, double a, double b){
    vector<int> histogram(bins);
    double step = (b-a) / bins;
    for(int i=0; i < signal.size(); i++){
        histogram[truncl((signal[i]-a)/step)]++;
        if(signal[i]==b){
            histogram[bins-1]++;
        }   
    }
    return histogram;
}

double mean_value(vector<int> histogram, int N_samples, int bins, double a, double b){
    double mean;
    double step = (b-a) / bins;
    for(int i=0; i< bins; i++){
        double temp = mean;
        mean = temp + (a+step*(2*i+1)/2)*((double)histogram[i]/N_samples);
    }
    cout << "average value: " << mean << "\n";
    return mean;
}

double variance(vector<int> histogram, double avg, int N_samples, int bins, double a, double b){
    double var;
    double step = (b-a) / bins;
    for(int i=0; i< bins; i++){
        double temp = var;
        var = temp + pow((a+step*(2*i+1)/2 - avg), 2)*((double)histogram[i]/N_samples);
    }
    cout << "variance: " << var << "\n";
    return var;
}

double quartile(vector<int> histogram, int N_samples, int bins, double a, double b, double num){
    int k=0;
    int tmp;
    double step = (b-a) / bins;
    for(int i = 0; i < bins; i++){
        tmp = k;
        k = tmp + histogram[i];
        if(k>=N_samples*num){
            return (a+step*(2*i+1)/2);
        }
    }
    return 0.;
}

vector<double> sine_signal(int N_samples, double T, double phi, double F=16.0, int quantization_levels_num=256, double quantization_min = -1., double quantization_max = 1.){
    vector<double> signal;
    for (int i = 0; i < N_samples; ++i){
        double signal_val = trunc((sin(2*3.1415*i/F/T + phi)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        signal.push_back(signal_val);
    }
    return signal;    
}

vector<double> triangular_signal(int N_samples, double T, double phi, double t, double F=16.0, int quantization_levels_num=256, double quantization_min = 0., double quantization_max = 1.){
    vector<double> signal;
    for (int i = 0; i < N_samples; ++i){
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
    for (int i = 0; i < N_samples; ++i){
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
        signal.push_back(*psv);
    }
    return signal;
}

bool is_power2(int x){
    for (int i = 1; i <= x; i *= 2){
        if (i == x) return 1;
    }
    return 0;
}

void preproc1(vector<double> &temp, int size, int sizeLog2){
    for(int i = size; i < pow(2, sizeLog2); i++){
        temp.push_back(0.);
    }    
}

void preproc2(vector<double> &temp, int size, int sizeLog2){
    double tmp;
    for(int i = size; i < pow(2, sizeLog2); i++){
        tmp=temp[i-size];
        temp.push_back(tmp);
    }    
}

void ComplexBitReverse(vector<complex<double>> &data, int size) {
    int middle = size / 2,
    revSize = size - 1,
    j = 0;
    for (int i = 0; i < revSize; ++i) {
        if(i < j) {
            swap(data[i], data[j]);
        }
        int k = middle;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

void FftDit( vector<complex<double>> &data,  int size, int sizeLog2, int dir ) {
    ComplexBitReverse(data, size);  // переставить в бит-реверсивном порядке
    int ptsInLeftDft,ptsInRightDft = 1;
    for ( int stage = 1; stage <= sizeLog2; ++stage )
    {
        ptsInLeftDft = ptsInRightDft;   // установить ptsInLeftDFT = 2**(stage-1)
        ptsInRightDft *= 2;             // установить ptsInRightDFT = 2**stage
        complex<double> twiddle = complex<double>(1.0, 0.0); // поворачивающий множ.
        double trigArg = M_PI / ptsInLeftDft;  
        // dir == 1 для прямого преобразования, dir == -1 для обратного
        complex<double> wFactor = complex<double>(cos(trigArg),-sin(trigArg)*dir);
        for( int butterflyPos = 0; butterflyPos < ptsInLeftDft; ++butterflyPos )
        {                             
            for(int topNode=butterflyPos; topNode < size; topNode+=ptsInRightDft )
            {                              
                int botNode = topNode + ptsInLeftDft;
                complex<double> temp = data[botNode] * twiddle;
                data[botNode] = data[topNode] - temp;
                data[topNode] += temp;
            }  // конец цикла по topNode

            twiddle *= wFactor;
        } // конец цикла "бабочка"
    } // конец цикла stage
}

int main() {
    double F = 16.;
    double N = 15;
    int N_samples = (int) N * F;
    double signal_value;
    double phi = 3.14 / 2.;
    double period = 1.;
    double t=0.5; //Длительность импульса
    double A=3.;
    int bins = 100;

    double a = -1.;
    double b = 1.;

    double s = 0.1;
    double P = 0.1;

    vector<complex<double>> data;
    vector<double> temp = sine_signal(N_samples, period, phi);
    //vector<double> temp = rectangular_signal(N_samples, period, phi, t, A);
    //vector<double> temp = triangular_signal(N_samples, period, phi, t);
    //vector<double> temp = random_signal(N_samples, t, A);
    if(!is_power2(N_samples)){
        //preproc1(temp, N_samples, trunc(log2(N_samples))+1);
        preproc2(temp, N_samples, trunc(log2(N_samples))+1);
        N_samples = temp.size();
        N = N_samples/F;
    }

    for(int i=0; i < temp.size(); i++){
        data.push_back(temp[i]);
    }
    FftDit(data, N_samples, (int) log2(N_samples), 1);
    for (int i = 0; i < N_samples / 2; i++) {
        swap(data[i], data[N_samples / 2 + i]);
    }
    vector_csv_complex(data, N, F);
    
    /*vector<double> signal = sine_signal(N_samples, period, phi);
    AWGN_signal(signal, s, a, b);
    Impulse_signal(signal, P, a, b); 
    vector_csv(signal, N, F);
    */
    
    /*vector<int> histogram = hist(signal, bins, a, b);

    for(int i = 0; i < histogram.size(); i++){
        cout << histogram[i] << " ";
    }
    double mean = mean_value(histogram, N_samples, bins, a, b);
    variance(histogram, mean, N_samples, bins, a, b);
    cout << "quartile distance: " << quartile(histogram, N_samples, bins, a, b, 0.75) - quartile(histogram, N_samples, bins, a, b, 0.25) << "\n";
    */
    return 0;
}