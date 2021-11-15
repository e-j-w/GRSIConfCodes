// g++ makeSimpleTable.cxx -o makeSimpleTable

#include <stdio.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

std::vector<std::string> chanNumber;
std::vector<std::string> hexColChan;
std::vector<std::string> channelName;

std::vector<std::string> collectors;
std::vector<std::string> digitizersHex;
std::vector<std::string> digitizers;
std::vector<std::string> digitizerChannelHex;

std::vector<std::string> digitizerChannel;

std::vector<std::string> digitizerADCName;


const std::string collector0[16] = {"grifadc129","grifadc131","grifadc132","grifadc133","grifadc134","grifadc135","grifadc140","grifadc141","grifadc108","grifadc109","grifadc106","grifadc111","grifadc112","grifadc113","grifadc114","grifadc115"};

const std::string collector1[16] = {"grifadc116","grifadc117","grifadc118","grifadc119","grifadc120","grifadc121","grifadc122","grifadc128","grifadc124","grifadc125","grifadc126","grifadc123","grifadc142","grifadc143","grifadc148","grifadc149"};

const std::string collector2[16] = {"grifadc150","grifadc151","grifadc152","grifadc153","null","grifadc137","grifadc138","grifadc139","null","grifadc145","grifadc154","grifadc147","grifadc100","grifadc101","grifadc102","grifadc104"};

const std::string collector3[16] = {"grifadc105","grifadc107","grifadc110","grifadc130","grifadc155","grifadc156","grifadc157","grifadc158","grifadc159","grifadc161","grifadc162","grifadc163","grifadc164","grifadc165","grifadc166","grifadc170"};

const std::string collector4[4] = {"grifadc173","grifadc178","grifadc136","grifadc177"};

const std::string collector5[1] = {"null"};

void readTable(const char* confFileName){
	
	ifstream inFile;
	inFile.open(confFileName);
	
	std::string inputColumn;
	
	if (inFile.is_open()){
		while (std::getline(inFile,inputColumn,'\t')){ //this starts a new line
			chanNumber.push_back(inputColumn); //save the first entry as the channel number
			std::getline(inFile,inputColumn,'\t'); //get the hex value
			hexColChan.push_back(inputColumn); //save the hex value
			std::getline(inFile,inputColumn,'\t'); //get the channel Name
			channelName.push_back(inputColumn); //save the channel name
			std::getline(inFile,inputColumn,'\n'); //get the rest of the line, which we don't care about
			/*std::getline(inFile,inputColumn,'\t'); //get data type (template) (which we don't care about)
			std::getline(inFile,inputColumn,'\t'); //get the gain (which we don't care about)
			std::getline(inFile,inputColumn,'\t'); //get the offset (which we don't care about)
			std::getline(inFile,inputColumn,'\t'); //get*/
		}
	}
	
	
	inFile.close();
	
	return;
}

void separateHex(std::vector<std::string> hexes){
	
	for (int i = 0; i < hexes.size(); i++){
		std::string val = hexes[i];
		std::string collector = val.substr(2,1);
		std::string digitizer = val.substr(3,1);
		std::string digChan  = val.substr(5,1);
		collectors.push_back(collector);
		digitizersHex.push_back(digitizer);
		digitizerChannelHex.push_back(digChan);
	}

	return;
}


void convertHexToDecimal(std::vector<std::string> hex, std::vector<std::string> &writeVector){ //shamelessly stolen from https://www.tutorialspoint.com/cplusplus-program-for-hexadecimal-to-decimal


	for (int i = 0; i < hex.size(); i++){
		std::string num = hex[i];
		int len = (num.size());
		int base = 1;
		int temp = 0;
		for (int j=len-1; j>=0; j--) {
		  if (num[j]>='0' && num[j]<='9') {
			 temp += (num[j] - 48)*base;
			 base = base * 16;
		  }
		  else if (num[j]>='a' && num[j]<='f') {
			 temp += (num[j] - 87)*base;
			 base = base*16;
		  }
		}
		//temp++; //increment temp by one because the ADC channels index at 1, not 0
		writeVector.push_back(std::to_string(temp));
	}

	return;
}

void determineDigitizer(std::vector<std::string> collectorNum, std::vector<std::string> digitizerNum){

	for (int i = 0; i < collectorNum.size(); i++){
		if ("0" == collectorNum[i]){
			digitizerADCName.push_back(collector0[std::stoi(digitizerNum[i])]);
		}
		else if ("1" == collectorNum[i]){
			digitizerADCName.push_back(collector1[std::stoi(digitizerNum[i])]);
		}
		else if ("2" == collectorNum[i]){
			digitizerADCName.push_back(collector2[std::stoi(digitizerNum[i])]);
		}
		else if ("3" == collectorNum[i]){
			digitizerADCName.push_back(collector3[std::stoi(digitizerNum[i])]);
		}
		else if ("4" == collectorNum[i]){
			digitizerADCName.push_back(collector4[std::stoi(digitizerNum[i])]);
		}
		else if ("5" == collectorNum[i]){
			digitizerADCName.push_back(collector5[std::stoi(digitizerNum[i])]);
		}
		else {
			digitizerADCName.push_back("-1");
		}
	}

	return;
}

void writeFile(std::vector<std::string> channels, std::vector<std::string> digitizerName, std::vector<std::string> digitizerChanNumb, std::vector<std::string> channelNaming){

	ofstream outFile;
	outFile.open("Conf_file_simple.txt");
	
	if (outFile.is_open()){
	
		outFile << "DAQ Chan\tADC Name\tADC Chan\tChan Name";
	
		for (int i = 0; i < channels.size(); i++){
			outFile << "\n" << channels[i] << "\t" << digitizerName[i] << "\t" << digitizerChanNumb[i] << "\t" << channelNaming[i];
		}
	}
	
	outFile.close();

	return;
	
}

int main(int argc, char * * argv){
	if (1 == argc){
		std::cout << "\nError: you need to give me a configuration file!";
		return 0;
	}
	if (2 == argc){
		readTable(argv[1]);
	}
	else{
		std::cout << "\nToo many inputs. Exiting now!";
		return 0;
	}
	
	separateHex(hexColChan);
	convertHexToDecimal(digitizerChannelHex,digitizerChannel);
	convertHexToDecimal(digitizersHex,digitizers);
	determineDigitizer(collectors,digitizers);
	
	writeFile(chanNumber,digitizerADCName,digitizerChannel,channelName);
	
	
	/*for (int i = 0; i < chanNumber.size(); i++){
		std::cout << "\nChannel number is " << chanNumber[i] << ", hex value is " << hexColChan[i] << ", channel name is " << channelName[i] << ", collector hex is " << collectors[i] << ", digitizer hex is " << digitizersHex[i] << ", digitizer channel hex is " << digitizerChannelHex[i] << ", digitizer channel is " << digitizerChannel[i] << ", digitizer name is " << digitizerADCName[i];
	}*/

	return 0;
}
