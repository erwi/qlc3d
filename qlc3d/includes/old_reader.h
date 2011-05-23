#ifndef READER_H
#define READER_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <locale>
#define READER_SUCCESS      1
#define READER_NOT_FOUND    -2   // variable not found
#define READER_BAD_VALUE    -3   //
#define READER_BAD_FORMAT   -4   // missing "=" etc.
#define READER_BAD_FILE	    -5   // file not open etc.
class Reader
{
    public:
		std::fstream file;


        Reader();
        virtual ~Reader();
        bool openFile(const std::string& filename);
        void closeFile();
	bool gotoLine(const unsigned int& l); // uggly seek line number in txt files

        int readString(std::string var, std::string& val);
	int readDlmLine(std::vector <double>& val, std::string dlm  );
	int findStringLine(std::string str,		// finds line consisting of string str,
			   unsigned int offset = 0);	// returns line number with respect to offset


	template <class T>
	int readNumberArray(std::string var, std::vector<T>& val)
	{

			using std::string;
			using std::cout;
			using std::endl;
			//cout << "readNumberArray:"<<var<<endl;

			val.clear();
			file.seekg(0);
			if (!file.good()) return READER_BAD_FILE;

			//file.clear(); // clear previous eof bits

            std::string line;

			while( ! file.eof() ){
					getline(file , line );

                    removeBlanks(line);
                    removeComments(line);
					//std::cout << "comparing "<< line << " with " << var << std::endl;
                    if (line.find(var) == 0){//std::string::npos){

                        size_t pos = line.find("=");    // find "=" sign
                        if (pos == std::string::npos)
                            return READER_BAD_FORMAT;

                        // get rid overything before and including "=" sign
                        line = line.substr(pos+1 , line.length() - (pos+1) );

                        // find opening and closing square braces "[" and "]"
                        size_t poso = line.find("[");
                        size_t posc = line.find("]");
                        if ( ( poso==string::npos ) || (posc == string::npos) || ( posc < ( poso+1)) )
                            return READER_BAD_FORMAT;

                        line = line.substr(poso+1 , posc - (poso+1) );

                        bool good = true;

                        while( good ){ // while valid numeric characters

                            pos = line.find_first_of(","); // first numeric character ends at pos

                            double temp_val;

                            if (std::stringstream( line.substr(0,pos) ) >> temp_val){ // test if numeric character
                                good = true;

                                val.push_back( (T) temp_val); // add value to vector with type casting

                                // Remove read characters from line
                                if (pos!= string::npos){// if not end of string
                                    line = line.substr(pos+1, line.length() - (pos+1) );
                                }
                                else
                                    return READER_SUCCESS; // exits when end of array is found

                            }// end if numeric character
                            else
                                good = false; // non-numeric character found
                            }// end while numeric characters
                        return READER_BAD_VALUE; // returned if illegal character found

                    }// end if variable found
            }//end while read file loop
			file.seekg(0);
			file.clear();
            return READER_NOT_FOUND;
	}// end readNumberArray

        template <class T>
int readNumber(std::string var, T& val){
	this->file.seekg(0);
	this->file.clear(); // clear eof bit
            std::string line;
            double temp_val; // temp value, will be typecast to int/float after reading
			while( ! file.eof() ){
				getline(file , line);
                //std::cout << "looking for " << var << " line is " << line << std::endl;
				removeBlanks(line);
                removeComments(line);
                //std::cout << "looking for " << var << " line is " << line << std::endl;
                if (line.find(var) == 0){// if var found in beginning of line
                    //std::cout << "var found !" << std::endl;

                    size_t pos = line.find ("=");
		    if (pos == std::string::npos){
			std::cout <<"bad format: "<< line << std::endl;
                        return READER_BAD_FORMAT;
		    }
			std::string s_val = line.substr(pos+1 , line.length() - (pos+1) );


			 if (std::stringstream(s_val) >> temp_val ){
				 val = (T) temp_val; // typecast to correct number format
				 return READER_SUCCESS;
			 }

		    else{
				std::cout <<"bad value: " << s_val << std::endl;
				return READER_BAD_VALUE;
			}
			}// if var found
            } // read file line by line loop
			file.seekg(0);
			file.clear();
            return READER_NOT_FOUND;
}
// end int readValue


        std::string getErrorString(int e);




    protected:
    private:
        bool removeBlanks(std::string& str);// removes spaces
	bool removeBlanksAtBeginning(std::string& str);// removes spaces at beginning of line only
	bool removeBlanksAtEnd(std::string& str);// removes empty spaces at end of string
	bool removeComments(std::string& str);// removes anythign after a "#" character
	void stringToLowercase(std::string& str); // converts string to all lowercase

};




#endif // READER_H
