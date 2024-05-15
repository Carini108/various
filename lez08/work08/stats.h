#include <iostream>
#include <fstream> 
//STL
#include <vector>  
//SORT SI STL
#include <algorithm>
//DEBUGGING
#include <assert.h>

using namespace std;

// ============================
// last update: may 14th 2024
// ============================
// FUNCTIONS 
// ============================


// VIDEO PRINT ===============================
template <typename T> void Print( const vector<T> & v ) {
  for (int i=0; i<v.size(); i++) cout << v[i] << endl;  
  return;
};

// FILE PRINT ================================
template <typename T> void Print( const vector<T> & v  , char* filename ) {
	ofstream fout(filename);
	if (!fout) {
			cerr << "Cannot create file " << filename<< endl;
			return;
	}
	for (int i=0; i<v.size(); i++) 
		if (i==0) { fout << v[i]; } else { fout <<"\n"<<v[i]; }
	fout.close();
	return;
};

// SMART READING ==========================================
// NEED NOT KNOW THE FILE LENGTH
template <typename T> vector<T> ReadAll( const char* filename ) {
	vector<T> v;
	ifstream fin(filename);
	if (!fin) {
		cout << "Cannot open file "<< filename << endl;
		exit(70);
	} else {
		while (!fin.eof()) {
			T val = 0;
			fin >> val;
			v.push_back(val);
		}
	}
	fin.close();
	return v;
};

// AVERAGE =========================================
template <typename T> double Average( const vector<T> & v) {
	double average = 0;
    if (v.size()==0) {
        return 0;
    } else {
        for ( int k = 0 ; k < v.size() ; k++ )
            average = static_cast<double>(k)/static_cast<double>(k+1)*average + v[k]/static_cast<double>(k+1); //sum of all data
        return average;
    }
};

// MEDIAN ==================================
template <typename T> double Median( vector<T> v ) {
    sort( v.begin(), v.end() ) ;  // Use the STL sort
	double median = 0;
    if (v.size()%2 == 0) //even vector size
        median = ( v[ v.size()/2-1 ] + v [ v.size()/2 ] )/2;
    else //odd vector size
        median = v[ v.size()/2 ];
    return median;
};

// VARIANCE ===========================================
template <typename T> double Variance( const vector<T> & v) {
    double average = Average (v);
    double variance = 0;
    if (v.size()==0) {
        return 0;
    } else {
        for ( int k = 0 ; k < v.size() ; k++ ){
            variance += (average-v[k])*(average-v[k]); //sum of all the squares
        }
        variance = variance / double (v.size());
        return variance;
    }
};
/*
// VARIANCE ===========================================
template <typename T> double Variance( const vector<T> & v) {
    double average = Average (v);
    double variance = 0;
    if (v.size()==0) {
        return 0;
    } else {
        double old_avg, avg = 0. ; 
        for ( int k = 0 ; k < v.size() ; k++ ){
            old_avg = avg;
            avg = static_cast<double>(k)/static_cast<double>(k+1)*average + 1./static_cast<double>(k+1)*v[k];
            variance = 1./static_cast<double>(k+1) * (static_cast<double>(k)*variance+v[k]*v[k]+static_cast<double>(k)*old_avg*old_avg)-avg*avg;
        }
        return variance;
    }*/