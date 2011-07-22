/*
 * parameterFile.h
 *
 *  Created on: Jul 22, 2011
 *      Author: fogelson
 */



#ifndef PARAMETERFILE_H_
#define PARAMETERFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

using namespace std;

namespace CFD{
	class ParameterFileReader{
		string hash;
		string filename;
		ifstream in;
		map<string, double> parameters;

	public:
		ParameterFileReader(){}
		ParameterFileReader(string filename){
			hash = "#";
			setFilename(filename);
			readFile();
		}
		void readFile(){
			in.open(filename.c_str(), ios::in);


			while(!in.eof()){
				string currentLine;
				getline(in,currentLine);

				int lineLength = currentLine.length();
				// Ignore anything after a # symbol
				for(int i = 0; i < lineLength; i++){
					if(currentLine.data()[i] == '#'){
						lineLength = i;
						break;
					}
				}
				// Ignore leading whitespace
				int leadingWhite = 0;
				for(int i = 0; i < lineLength; i++){
					if(currentLine.data()[i] == ' ' || currentLine.data()[i] == '\t'){
						leadingWhite++;
					}
					else{
						break;
					}
				}

				stringstream paramName;

				// Try to read parameter name
				for(int i = leadingWhite; i < lineLength; i++){

				}
			}

			in.close();
		}
		~ParameterFileReader(){
			if(in.is_open()){
				in.close();
			}
		}
		string getFilename(){
			return filename;
		}
		void setFilename(string filename){
			this->filename = filename;
		}
	};
}
#endif /* PARAMETERFILE_H_ */
