// -*- c++ -*-
///////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2001 Oh-Wook Kwon, all rights reserved. ohwook@yahoo.com
//
//                          Easy Matrix Template Library
// 
// This Easy Matrix Template Library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the authors shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
// 
// Everyone is granted permission to copy, modify and redistribute this
// Easy Matrix Template Library, provided:
//  1.  All copies contain this copyright notice.
//  2.  All modified copies shall carry a notice stating who
//      made the last modification and the date of such modification.
//  3.  No charge is made for this software or works derived from it.  
//      This clause shall not be construed as constraining other software
//      distributed on the same medium as this software, nor is a
//      distribution fee considered a charge.
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Filename: MatrixUtil.h
///////////////////////////////////////////////////////////////////////////////
#ifndef _MATRIX_UTIL_H_
#define _MATRIX_UTIL_H_

#include "./Matrix.h"
#include "./Random.h"

template <typename T> class Vector;
template <typename T> class Matrix;


inline string Int2str(int i) { 
	string s; char x[128];
	sprintf(x,"%d",i); s=x;
	return s;
}


inline string Num2str(int i) {
	string s; char x[128];
	sprintf(x,mtl_format_int,i); s=x;
	return s;
}


inline string Num2str(double f) {
	string s; char x[128];
	sprintf(x,mtl_format_float,f); s=x;
	return s;
}


template <typename T> string Num2str(const complex<T>& c) {
	char x[128], re[64], im[64], sum[132];
	int i;
	sprintf(re,mtl_format_float,c.real());
	sprintf(im,mtl_format_float,c.imag());
	strcpy(sum,"(");
	for(i=0;i<strlen(re);i++) if(re[i] != ' ') break;
	strcat(sum,&re[i]);
	strcat(sum,",");
	for(i=0;i<strlen(im);i++) if(im[i] != ' ') break;
	strcat(sum,&im[i]);
	strcat(sum,")");
	sprintf(x,mtl_format_complex,sum);
	string s=x;
	return s;
}


template <typename T> int Print(ostream& s, const char *name, T val){
	if(name != 0 && name[0] != 0) s << name << " ";
	s << 0 << " " << val << endl;
	return 1;
}


inline int Print(ostream& s, const string& name, int val){
	return Print(s,name.c_str(),val);
}


inline int Print(ostream& s, const string& name, double val){
	return Print(s,name.c_str(),val);
}


template <typename T> int Read(istream& s, const char *name, T& val){
	string dummy;
	s >> dummy;
	if(name != 0 && name[0] != 0) {
		assert(name==dummy);
	}
	int itmp;
	s >> itmp;
	assert(itmp==0);
	s >> val;
	return 1;
}


inline int Read(istream& s, const string& name, int& val){
	return Read(s, name.c_str(), val);
}


inline int Read(istream& s, const string& name, double& val){
	return Read(s, name.c_str(), val);
}


inline Vector<string> mtl_split(const string& items, const char* delimiter);


inline Vector<string> mtl_split(const string& items) {
	return mtl_split(items, " \t"); // to avoid internal compiler error in Visual C++ 5.0
}


inline bool mtl_isspace(int c) {
	if(c == ' ' || c == '\t' || c == '\r' || c == '\n') return true;
	else return false;
}


inline string& mtl_trimspace(string& word)
{
	int i;
	int n=word.length();
	for(i=n-1; i>=0; i--){
		if(!mtl_isspace(word[i])) break;
	}
	if(n-i-1>0) word.erase(i+1,n-i-1);
	n=i+1;
	for(i=0;i<n;i++){
		if(!mtl_isspace(word[i])) break;
	}
	if(i>0)	word.erase(0,i);
	return word;
}


inline Vector<string> mtl_split(const string& items, const char* delimiter) {
	Vector<string> itemArray;
	string token("");
	int i;	
	for(i=0; i<items.size(); ++i) {
		if(strchr(delimiter,items[i]) != NULL) {
			if(token.size() == 0) continue;
			else {
				itemArray.Add(token);
				token = "";
			}
		}
		else token += items[i];
	}
	if(token != "") itemArray.Add(token);
	for(i=1; i<=itemArray.Size(); ++i) {
		mtl_trimspace(itemArray[i]);
	}
	return itemArray;	
}


#endif
