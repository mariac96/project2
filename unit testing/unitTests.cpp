#include "acutest.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <bits/stdc++.h>
#include <numeric>
#include <map>
#include "main.hpp"
#include "lsh.hpp"
#include "frechet2.hpp"
#include <math.h>
using namespace std;


void test_binaryToDecimal(void){

  string temp = "010100001111110101100";

  int result = binaryToDecimal(temp);
  TEST_CHECK(result == 663468);
  TEST_MSG("Result1: %d", result);

  string temp2 = "010101100110";

  int result2 = binaryToDecimal(temp2);
  TEST_CHECK(result2 == 1382);
  TEST_MSG("Result2: %d", result2);
}

void test_compareDouble(void){
  double a = 42.5;
  double b = 42.50000001;
  TEST_CHECK(compareDouble(a,a));
  TEST_CHECK(!compareDouble(a,b));

}

void test_euclidModulo(){
  TEST_CHECK(euclidModulo(5,3)==2);
  TEST_CHECK(euclidModulo(88,59)==29);
  TEST_CHECK(euclidModulo(68,-20)==-12);
}


void test_euclidModuloLong(){
  TEST_CHECK(euclidModuloLong(5,309868)==5);
  TEST_CHECK(euclidModuloLong(88,5912647)==88);
  TEST_CHECK(euclidModuloLong(-68,-201240184)==-68);
}

TEST_LIST = {
    { "test_binaryToDecimal", test_binaryToDecimal },
    { "test_compareDouble", test_compareDouble },
    { "test_euclidModulo", test_euclidModulo },
    { "test_euclidModuloLong", test_euclidModuloLong },
    { NULL, NULL }
};
