#include "Curvature-new.h"
#include "MEMORY.H"
#include <cassert>
#include <iostream>

using namespace std;

const int ARRAYSIZE = 125;
const int ARRAYLEN = 5;

void run_tests();
void assign_iPos(DMATRIX d);
void assign_jPos(DMATRIX d);
void assign_kPos(DMATRIX d);
void assign_fracSolid(DMATRIX d);
void assign_slushdir(DMATRIX d);

int main(){
    DMATRIX ipos = (DMATRIX)Create3DArray(5,5,5,sizeof(double));
    DMATRIX jpos= (DMATRIX)Create3DArray(5,5,5,sizeof(double)); 
    DMATRIX kpos= (DMATRIX)Create3DArray(5,5,5,sizeof(double)); 
    DMATRIX fracsolid= (DMATRIX)Create3DArray(5,5,5,sizeof(double)); 
    DMATRIX slushdir= (DMATRIX)Create3DArray(5,5,5,sizeof(double)); 
    assign_iPos(ipos);
    assign_jPos(jpos);
    assign_kPos(kpos);
    assign_fracSolid(fracsolid);
    assign_slushdir(slushdir);
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                assert(ipos[i][j][k] == static_cast<double>(i+j+k));
                assert(jpos[i][j][k] == static_cast<double>(i-j+k));
                assert(kpos[i][j][k] == static_cast<double>(i*j+k));
                assert(fracsolid[i][j][k] == static_cast<double>(j+k));
                assert(slushdir[i][j][k] == static_cast<double>(j+k*i));
            }
        }
    }
	char c = 'a';
	while(c != 'q'){
		cin >> c;
		cout << "No errrors!" << endl;
	}
    return 0;
}

void assign_iPos(DMATRIX d){
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                d[i][j][k] = static_cast<double>(i+j+k);
            }
        }
    }
}

void assign_jPos(DMATRIX d){
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                d[i][j][k] = static_cast<double>(i-j+k);
            }
        }
    }
}

void assign_kPos(DMATRIX d){
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                d[i][j][k] = static_cast<double>(i*j+k);
            }
        }
    }
}

void assign_fracSolid(DMATRIX d){
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                d[i][j][k] = static_cast<double>(j+k);
            }
        }
    }
}

void assign_slushdir(DMATRIX d){
    for(int i = 0; i < ARRAYLEN; i++){
        for(int j = 0; j < ARRAYLEN; j++){
            for(int k = 0; k < ARRAYLEN; k++){
                d[i][j][k] = static_cast<double>(j+k*i);
            }
        }
    }
}
