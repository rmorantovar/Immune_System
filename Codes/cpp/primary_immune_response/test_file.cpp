//
//  test_file.cpp
//  
//
//  Created by Roberto Moran Tovar on 27.05.21.
//
#include "../lib/Immuno_functions.hpp"

#include <stdio.h>


int main(int argc, char* argv[]) //argv has ;
{
    ofstream fout ("../Input_files/MM.txt");
    
    for (int j=0; j<20; j++) {
        for (int i = 0; i<20; i++) {
            
            if (i==j) {
                fout << -1 << "\t";
            } else {
                fout << 0 << "\t";
            }
        }
        fout << endl;
    }
    
}
