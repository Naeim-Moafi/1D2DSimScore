#include <cmath>
#include <cstdio>
#include <exception>
#include "SimilarityScores.h"
#include <stdexcept>
#include <iostream>


//initializing ctor
SimilarityScores::SimilarityScores(int TP, int TN, int FP, int FN)
: TP(TP), TN(TN), FP(FP), FN(FN)
{}


// setCategories
void SimilarityScores::setCategories(int TP, int TN, int FP, int FN)
{
	this->TP = TP;
	this->TN = TN;
	this->FP = FP;
	this->FN = FN;
}

// calcMCC
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcMCC() const
{
	double numerator = ((static_cast<double>(TP) * static_cast<double>(TN)) - (static_cast<double>(FP) * static_cast<double>(FN)));
	double denominator = sqrt((static_cast<double>(TP + FP) * static_cast<double>(TP + FN) * static_cast<double>(TN + FP) * static_cast<double>(TN + FN)));
	
	double mcc = numerator / denominator;
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return mcc;
}

// calcPrecision
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcPrecision() const
{
	double numerator = static_cast<double>(TP);
	double denominator = static_cast<double>(TP + FP);
	double ppv = numerator /denominator;
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return ppv;
	
}

// calcRecall
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcRecall() const
{
	double numerator = static_cast<double>(TP);
	double denominator = static_cast<double>(TP + FN);
	double recall = numerator / denominator;
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return recall;
}

// calcJIndex
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcJIndex() const
{
	double numerator = static_cast<double>(TP);
	double denominator = TP + FP + FN;
	double jIndex = numerator / denominator;
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return jIndex;
}

// calcFscore
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcFscore() const
{
	double numerator = calcPrecision() * calcRecall();
	double denominator = calcPrecision() + calcRecall();
	double fScore = 2 * numerator / denominator;
	
	
	if(denominator == 0 || calcPrecision() == -1.1 || calcRecall() == -1.1)
	{
		return -1.1;
	}
	
	return fScore;
}

// calcFMIndex
// using nodiscard attribute to make sure the return value is not being ignored. 
double SimilarityScores::calcFMIndex() const
{
	if(calcPrecision() == -1.1 || calcRecall() == -1.1) return -1.1;
	return sqrt(calcPrecision() * calcRecall());
}

		
// Specificity
double SimilarityScores:: calcSpecificty() const
{
	double numerator = static_cast<double>(TN);
	double denominator = static_cast<double>(TN) + static_cast<double>(FP);
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return (numerator/denominator);
}

// Balance Accuracy
double SimilarityScores:: calcBA() const
{
	if(calcSpecificty() == -1.1 || calcRecall() == -1.1) return -1.1;
	return (calcSpecificty() + calcRecall()) * 0.5;
}

// False omision rate
double SimilarityScores:: calcFOR() const
{
	double numerator = static_cast<double>(FN);
	double denominator = static_cast<double>(FN) + static_cast<double>(TN); 
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return (numerator / denominator);
}

// Prevalence threshold
double SimilarityScores:: calcPT() const
{
	double numerator = 1 - calcRecall();
	double denominator = sqrt(numerator) + sqrt(1 - calcSpecificty());
	
	if(denominator == 0 || calcRecall() == -1.1)
	{
		return -1.1;
	}
	
	return (numerator / denominator);
}

// Critical success index
double SimilarityScores:: calcCSI() const
{
	double numerator = static_cast<double>(TP);
	double denominator = static_cast<double>(TP) + static_cast<double>(FN) + static_cast<double>(FP);
	
	if(denominator == 0)
	{
		return -1.1;
	}
	
	return (numerator / denominator);
}

// Markedness
double SimilarityScores:: calcMK() const
{
	double NPV = 1 - calcFOR();
	if(calcPrecision() == -1.1 || calcFOR() == -1.1) return -1.1;
	return (calcPrecision() + NPV - 1);
}



//#define _MAIN_
#ifdef _MAIN_

using namespace std;

int main()
{
	SimilarityScores sim1(1, 2, 3, 4);
	SimilarityScores sim2(4, 3, 2, 1);
	
	cout << "caculate mcc for sim1: " <<  sim1.calcMCC() << endl;
	cout << "caculate mcc for sim2: " <<  sim2.calcMCC() << endl;
	
	cout << "caculate ppv for sim1: " <<  sim1.calcPrecision() << endl;
	cout << "caculate ppv for sim2: " <<  sim2.calcPrecision() << endl;
	
	cout << "caculate recall for sim1: " <<  sim1.calcRecall() << endl;
	cout << "caculate recall for sim2: " <<  sim2.calcRecall() << endl;
		
	cout << "caculate Jindex for sim1: " <<  sim1.calcJIndex() << endl;
	cout << "caculate Jindex for sim2: " <<  sim2.calcJIndex() << endl;
	
	cout << "caculate Fscore for sim1: " <<  sim1.calcFscore() << endl;
	cout << "caculate Fscore for sim2: " <<  sim2.calcFscore() << endl;
	
	cout << "calculate FMIndex for sim1:" << sim1.calcFMIndex() << endl;
	cout << "calculate FMIndex for sim2:" << sim2.calcFMIndex() << endl;
	
	sim1.setCategories(3, 4, 1, 2);
	sim2.setCategories(4, 1, 3, 2);
	
	cout << "\ncaculate mcc for sim1: " <<  sim1.calcMCC() << endl;
	cout << "caculate mcc for sim2: " <<  sim2.calcMCC() << endl;
	
	cout << "caculate ppv for sim1: " <<  sim1.calcPrecision() << endl;
	cout << "caculate ppv for sim2: " <<  sim2.calcPrecision() << endl;
	
	cout << "caculate recall for sim1: " <<  sim1.calcRecall() << endl;
	cout << "caculate recall for sim2: " <<  sim2.calcRecall() << endl;
		
	cout << "caculate Jindex for sim1: " <<  sim1.calcJIndex() << endl;
	cout << "caculate Jindex for sim2: " <<  sim2.calcJIndex() << endl;
	
	cout << "caculate Fscore for sim1: " <<  sim1.calcFscore() << endl;
	cout << "caculate Fscore for sim2: " <<  sim2.calcFscore() << endl;
	
	cout << "calculate FMIndex for sim1:" << sim1.calcFMIndex() << endl;
	cout << "calculate FMIndex for sim2:" << sim2.calcFMIndex() << endl;
	
	return 0;
}


#endif //_MAIN_
