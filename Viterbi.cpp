# include <iostream>
# include <string>
# include <list>
# include <vector>

using namespace std;
//This function is to get the best frontier nodes.
//Input:	regionLoss: the loglikelihood for all regionsm
//			elevation: elevation information
//Output:	result: the vector of the frontier nodes
//Global:	lambda: a hyperparameter for hight difference(i.e. a parameter of second term in loss function formula).
vector<size_t> Viterbi(vector<vector<double>>* regionLoss, vector<vector<double>>* elevation)
{
	vector<double> regionStep1;
	vector<double> regionStep2;
	vector<double> elevation1;
	vector<double> elevation2;
	vector<vector<double>> delta;
	vector<vector<size_t>> path;
	vector<double>* tempLoss;
	double loss, tempBest,max;
	size_t individualPath, sta;
	vector<size_t> result;
	for (size_t i =0;i<regionLoss->size()-1;i++)
	{
		regionStep1 = regionLoss->at(i);
		regionStep2 = regionLoss->at(i+1);
		elevation1 = elevation->at(i);
		elevation2 = elevation->at(i+1);
		for (size_t j =0;j<regionStep2->size();j++)
		{
			vector<double> temp;
			vector<size_t> tempPath;
			tempBest = 0;
			
			for (size_t k =0;k<regionStep1->size();k++)
			{
				if (i == 0)
					loss = regionStep1.at(k) + parameter.lambda * abs(elevation1.at(i) - elevation2.at(i + 1));
				else
				{
					tempLoss = &delta.back();
					loss = tempLoss->at(k) + regionStep1.at(k) + parameter.lambda * abs(elevation1.at(i) - elevation2.at(i + 1));
				}
				if (loss > tempBest)
				{
					tempBest = loss;
					individualPath = k; //Save the best route and remove all other routes
				}
			}
			temp.push_back(tempBest);
			tempPath.push_back(individualPath); //Save the best route in one region
			
		}
		delta.push_back(temp);
		path.push_back(tempPath);
	}
	//Find the path
	max = delta.back().at(0);
	sta = 0;
	for (size_t i = 0; i < delta.back().size(); i++)
	{
		if (delta.back().at(i) > max)
		{
			max = delta.back().at(i);
			sta = i;
		}
	}
	result.push_back(sta);
	for (size_t i = 0; i < regionLoss->size(); i++)
	{
		sta = path.at(regionLoss->size() - i - 1).at(sta);
		result.push_back(sta);
	}
	reverse(result.begin(), result.end());
	return result;
}
