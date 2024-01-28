#include <iostream>
#include <fstream>
#include <string> 
#include "matrix.cpp"
#include "pbPlots.hpp"
#include "supportLib.hpp"
using namespace std;

int main(){
    int N;
    double L, K, Dt, T0, Tn, k;
    L = 2;
    N = 100;
    K = 200;
    Dt = 0.01;
    T0 = 50;
    Tn = 500;
    k = -0.08;

    double Dx, F;
    Dx = L/N;
    F = (k*Dt)/(Dx*Dx);

    Matrix A(N-1,N-1);
    for(int i=0; i<N-1; i++){
        for(int j=0; j<N-1; j++){
            if(i == j){
                A(i,j) = 1 - (2*F);
                if(i == 0 && j == 0){
                    A(i,j+1) = F;
                }else if(i == N-2 && j == N-2){
                    A(i,j-1) = F;
                }else{
                    A(i,j-1) = F;
                    A(i,j+1) = F;
                }
            }else{
                if(A(i,j) != F){
                    A(i,j) = 0;
                }
            }
        }
    }

    Matrix T(N+1,K+1);         //mudei
    for(int i=0; i<K+1; i++){  //mudei
        T(0,i) = T0;
        T(N,i) = Tn;
    }
    for(int i=1; i<N; i++){
        T(i,0) = 0;  
    }

   // FILE* file = fopen("alicia.datas", "w");
    
    Matrix B(N-1, 1);
    Matrix X(N-1, 1);
    int aux = -1;
    for(double i=0; i<K+1; i+=Dt){  //mudei
        if(i >= (aux+1)){
            aux++;
            for(int j=1; j<N; j++){
                X(j-1,0) = T(j,aux);
            }
        }
        X(0,0) += T0;
        X(N-2,0) += Tn;
     /*   for(double j=1; j<N; j++){
            if(j == 1){
                B(j-1,0) = T(j,i) + T0;
            }else if(j == N-1){
                B(j-1,0) = T(j,i) + Tn;
            }else{
                B(j-1,0) = T(j,i);
            }
        }*/
        //X = A.inverse()*B;
        X = A.inverse()*X;
        for(int j=1; j<N; j++){
            T(j,aux+1) = X(j-1,0);  //mudei i -> aux
        }
    }
    cout << T;
    cout << endl;
    /*double a, b, c, d, e, f, g, h, i, j;
    double sum = 0;
    for(int k=0; k<10; k++){
        for(int l=0; l<N; l++){
            sum += T(l, k);
        }
        sum = sum/N+1;
        switch (k){
        case 0:
            a = sum;
            cout << a << endl;
            break;
        case 1:
            b = sum;
            cout << b << endl;
            break;
        case 2:
            c = sum;
            cout << c << endl;
            break;
        case 3:
            d = sum;
            cout << d << endl;
            break;
        case 4:
            e = sum;
            cout << e << endl;
            break;
        case 5:
            f = sum;
            cout << f << endl;
            break;
        case 6:
            g = sum;
            cout << g << endl;
            break;
        case 7:
            h = sum;
            cout << h << endl;
            break;
        case 8:
            i = sum;
            cout << i << endl;
            break;
        case 9:
            j = sum;
            cout << j << endl;
            break;
        default:
            break;
        }
        sum = 0;
    }*/

///////////////////////////////////////////////////////////////////////////GRÁFICO////////////////////////////////////////////////////////////////////////
    bool success;
	StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
	RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();

    double ret[N+1];            //vetor que recebe todos os valores da matriz
    vector<vector<double>> ys;  //vetor de vetores da classe vector (vai armazenar os valores da matriz dividido pelo tempo K)
    for(int j = 0; j <=K; j++){
        for(int i = 0; i<=N; i++){
            if(j%20 == 0){
                ret[i] = T(i,j);
            }
           // cout << ret[i] << " ";
        }
        vector<double> retCopy(ret, ret+sizeof(ret)/sizeof(double));   //vetor da classe vector que irá copiar double ret
        ys.push_back(retCopy);  //valores das linhas da matriz Θ^k no eixo Y
       // cout << endl;
    }
    vector<double> xs;      //eixo X recebe apenas os valores de N (total de partes da barra)
    for(int n=0; n<N+1; n++){
        xs.push_back(n);
    }

    vector<ScatterPlotSeries*> series;      //esse vetor armazena todas as series que deverão ser plotadas no grafico
    //essa quantidade é proporcional a divisao do tempo K

    for(int i=0; i<N+1; i++){
        ScatterPlotSeries *plot = GetDefaultScatterPlotSeriesSettings();
        plot->xs = &xs;
        plot->ys = &ys[i];              //variação dos valores em θ^K
        plot->linearInterpolation = true;
        plot->lineType = toVector(L"solid");
        plot->lineThickness = 2;
        plot->color = CreateRGBColor(0, 1, 0);
        series.push_back(plot);
    }
	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 600;
	settings->height = 400;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
    for(int i=0; i<K; i++){     //divisao de tempo K define a quantidade de series no grafico
        settings->scatterPlotSeries->push_back(series[i]);
    }

	success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
	if(success){
		vector<double> *pngdata = ConvertToPNG(imageReference->image);
		WriteToFile(pngdata, "datas.png");
		DeleteImage(imageReference->image);
	}else{
		cerr << "Error: ";
		for(wchar_t c : *errorMessage->string){
			cerr << c;
		}
		cerr << endl;
	}
	FreeAllocations();
	return success ? 0 : 1;
}

