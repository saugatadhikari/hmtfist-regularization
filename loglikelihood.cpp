#define _USE_MATH_DEFINES
#define LOGZERO -INFINITY
#define Dim 3 //input data dimension, 3 means no texture, 6 means texture included
#define cNum 2


#include<queue>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<numeric>
#include<vector>
#include<string>
#include <chrono>
#include<ctime>
#include<cmath>
#include<limits>
#include <cstdio>
//Here are extended functions for ln and e^().
double eexp(double x) {
	if (x == LOGZERO) {
		return 0;
	}
	else {
		return exp(x);
	}
}
double eln(double x) {
	if (x == 0) {
		return LOGZERO;
	}
	else if (x > 0) {
		return log(x);
	}
	else {
		cout << "Negative input error " << x << endl;
		int tmp;
		cin >> tmp;
	}
}
//As we discused, the parents cannot affect the probability of the frontier node(i.e. all parent nodes are class 1). I removed the detection of the number of parents.
// Input:	waterlikelihood: the probabilty for water node from Unet.(Not loglikelihood) The order should follow the tree sturcture.
//			drylikelihood: the probabilty for water node from Unet.The order should follow the tree sturcture.
//			treeInput: the Id of the node in tree structure. The are only nodes in main branch. The order of the nodes is from lower to higher(i.e. the first Id should be the id of the lowest node in one region, the river node).
//			Nodes: The information for parents and children.
//			treeLength: the length of the tree structure for traversing the tree structure. 
// Output:	The vector of the loglikelihood for every frontier nodes in one region.
// Global parameter: rho and Pi. In code, the names are parameter.Pi and parameter.rho. Please modify the names of these parameters
vector<double> getLoglikelihood(vector<double>* waterlikelihood, vector<double>* drylikelihood, sTree treeInput, vector<Node*> Nodes, size_t treeLength)



{

	double initialLog = 0;

	//Calculate Gain
	double curWaterProb, curDryProb,  curMaxGain = 0;
	vector<double> loglikelihood;

	double curGain = 0;

	//define temporary variables
	int curIdx = -1, newIdx = -1;
	Node* curNode = NULL;


	//get initial loglikelihood. All nodes are class 0.
	for (size_t i = 0; i < treeLength; i++) {//Traverse all nodes in tree structure from bottom to top
		curIdx = treeInput.nodeLevelIndexPair[i].second; //curIdx is the id of current node
		curNode = Nodes.at(curIdx);

		if (data.missingIndicator.at(newIdx) == 0) { // missing data. Give default probability
			curDryProb = eln(0.5);
			curWaterProb = eln(0.5);
		}
		else {
			curDryProb = eln(drylikelihood->at(curIdx)); // log form
			curWaterProb = eln(waterlikelihood->at(curIdx));
		}


		curNode->visited = true;
		parentsAllVisited = true;

		if (curNode->parents.size() == 0) { //node in river. The second term in formula is missing, so the second term is removed.

			if (curNode->next != NULL) {
				curGain = curDryProb + eln(1 - parameter.Pi);
			}
			else {//false node which is seperated from other nodes, chain length = 1, may due to mosaic or data error
				curGain = 0;
			}

		}
		else { // internal node
			if (curNode->next != NULL) {
				curGain = curDryProb + eln(1 - parameter.Pi);
			}
			else { //The highest node
				curGain = curDryProb + eln(1 - parameter.Pi);
			}
		}

		initialLog += curGain + initialLog;

	}

	curMaxGain = 0;


	curIdx = -1;
	newIdx = -1;
	curNode = NULL;
	//get other loglikelihood for different frontier nodes.
	for (size_t i = 0; i < treeLength; i++) {
		curIdx = tree.nodeLevelIndexPair[i].second;
		curNode = data.allNodes.at(curIdx);
		newIdx = data.old2newIdx.at(curIdx);

		if (data.missingIndicator.at(newIdx) == 0) { // missing data
			curDryProb = eln(0.5);
			curWaterProb = eln(0.5);
		}
		else {
			curDryProb = eln(drylikelihood->at(curIdx)); // log form
			curWaterProb = eln(waterlikelihood->at(curIdx));
		}


		curNode->visited = true;
		parentsAllVisited = true;

		//if (curWaterProb + log(1 - parameter.M) >= curDryProb + log(parameter.M)) {
		//	curWaterProb += log(1 - parameter.M);
		//}
		//else {
		//	curWaterProb = curDryProb + log(parameter.M);
		//}

		if (curNode->parents.size() == 0) { 

			if (curNode->next != NULL) {
				curGain = curWaterProb - curDryProb - eln(parameter.Pi) + eln(1 - parameter.Pi);
			}
			else {//false node which is seperated from other nodes, chain length = 1, may due to mosaic or data error
				curGain = 0;
			}

		}
		else { // internal node
			if (curNode->next != NULL) {
				curGain = curWaterProb - curDryProb + eln(parameter.rho) - eln(1 - parameter.rho) - eln(parameter.Pi) + eln(1 - parameter.Pi);
			}
			else { 
				curGain = curWaterProb - curDryProb + eln(parameter.rho) - eln(1 - parameter.rho) - eln(parameter.Pi) + eln(1 - parameter.Pi);
			}
		}


		if (curNode->parents.size() == 0) { // bottom node

			curMaxGain = initialLog += curGain;
			loglikelihood.push_back(curMaxGain);

		}
		else {
			curMaxGain = curMaxGain + curGain;
			loglikelihood.push_back(curMaxGain);
		}

	}
	return loglikelihood;

};