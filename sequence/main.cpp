#include <iostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <fstream>
#include <set>
#include <sstream>
#include "../server/server.hpp"
#include <algorithm>



using std::pair;
using std::string;
using std::vector;
using std::cout;

const string SEQUENCE_SERVICE_PREFIX = "/range?";
const string QUERY_DELIM = ";";
const string FIELD_DELIM = ",";




vector<string> splitString(const string& sourceString, const string& delimiter, bool dropEmptyChunks = true) {
	vector<string> out;
	int startPos = 0;
	while (true) {
		int delimPos = sourceString.find(delimiter, startPos);
		if (delimPos == string::npos) {
			break;
		}
		if (delimPos != startPos || !dropEmptyChunks) {
			out.push_back(sourceString.substr(startPos, delimPos - startPos));
		}
		startPos = delimPos + 1;
	}
	if (startPos != sourceString.size() || !dropEmptyChunks) {
		out.push_back(sourceString.substr(startPos, sourceString.size() - startPos));
	}
	return out;
}


char chrom2chromCode(const string& chromosome) {
	const char X_CHROM_CODE = 255;
	const char Y_CHROM_CODE = 254;
	int chromCode = 0;
	char buffer[3];
	{
		if (chromosome.c_str()[0] == 'x' || chromosome.c_str()[0] == 'X') {
			chromCode = X_CHROM_CODE;
		} else if (chromosome.c_str()[0] == 'y' || chromosome.c_str()[0] == 'Y') {
			chromCode = Y_CHROM_CODE;
		} else {
			buffer[0] = chromosome.c_str()[0];
			buffer[1] = chromosome.c_str()[1];
			buffer[2] = 0;
			chromCode = atoi((const char*)buffer);
		}
	}
	return chromCode;
}

struct TChromInterval {
	int ChromCode;
	int StartBP, EndBP;
	TChromInterval(const string& chromosome, const int start, const int end) :
		ChromCode(chrom2chromCode(chromosome)),
		StartBP(start),
		EndBP(end) {
	}
};




vector<TChromInterval> parseQuery(const http::server::request& req) {
	vector<TChromInterval> queries;
	{//parse
		string query = req.uri.substr(SEQUENCE_SERVICE_PREFIX.size(), req.uri.size() - SEQUENCE_SERVICE_PREFIX.size());
		vector<string> chunks = splitString(query, QUERY_DELIM);
		for (vector<string>::const_iterator chunkIt = chunks.begin(); chunkIt != chunks.end(); ++chunkIt) {
			vector<string> queryFields = splitString(*chunkIt, FIELD_DELIM);
			try {
				string queryChromosome = queryFields[0];
				int startBP = boost::lexical_cast<int>(queryFields[1]);
				int endBP = boost::lexical_cast<int>(queryFields[2]);
				queries.push_back(TChromInterval(queryChromosome, startBP, endBP));
			} catch (std::exception& e) {
				std::cerr << "exception: " << e.what() << "\n";
			}
		}
	}
	return queries;
}


const string DATA_FOLDER = "../data/sequence/chroms/";
const string DATA_FILE_FIELD_DELIM = " ";


struct TStructCode {
	bool Code[3];
};


const char BASE_CODES[5] = {'N', 'A', 'C',  'G', 'T'  };

char code2Base(const bool* code) {
	char index = (code[0] << 2) + (code[1] << 1) + code[2];
	return BASE_CODES[index];
}

char base2Code(char base) {
	switch(base) {
		case 'N': return 0;
		case 'n': return 0;

		case 'A': return 1;
		case 'a': return 1;

		case 'C': return 2;
		case 'c': return 2;

		case 'G': return 3;
		case 'g': return 3;

		case 'T': return 4;
		case 't': return 4;

		default: return 0;
	}

}



const int WELL_PACKED_VECTOR_SIZE = 16777216;
bool DEBUG = false;


class TSequenceDB {
public:
	TSequenceDB() {
		cout << "started\n";
		vector<string> chromosomes;
		for (int chrId = 1; chrId < 23; ++chrId) {
			chromosomes.push_back(boost::lexical_cast<std::string>(chrId));
		}
		chromosomes.push_back("X");
		chromosomes.push_back("Y");
		for (int chrIndex = 0; chrIndex < chromosomes.size(); ++chrIndex) {
			string dataFile = DATA_FOLDER + string("chr") + chromosomes[chrIndex] + ".fa";
			cout << dataFile << "\n";
			std::ifstream file(dataFile.c_str());
			if (!file.is_open()) {
				continue;
			}
			cout << dataFile << "\n";

			int vectorPosition = 0;
			{
				TChromInterval interval(chromosomes[chrIndex], 0, 0);
				Indices.push_back(pair<TChromInterval, vector<bool> >(interval, vector<bool>(WELL_PACKED_VECTOR_SIZE, false)));
			}
			vector<bool>* currentIndexPtr = &Indices.rbegin()->second;
			string line;
			int processed = 0;
			while (getline(file, line)) {
				{
					if (processed % 1000000 == 0) {
						cout << "..processed: " << processed << "\n";
					}
					++processed;
				}
				if (line.c_str()[0] == '>') {
					continue;
				}

				int toAdd = line.length() * 3;
				if (vectorPosition + toAdd >= WELL_PACKED_VECTOR_SIZE) {
					{
						int basesInChunk = vectorPosition / 3;
						int endBaseIndex =  Indices.rbegin()->first.StartBP + basesInChunk;
						Indices.rbegin()->first.EndBP = endBaseIndex;
						TChromInterval interval(chromosomes[chrIndex], endBaseIndex, 0);
						Indices.push_back(pair<TChromInterval, vector<bool> >(interval, vector<bool>(WELL_PACKED_VECTOR_SIZE, false)));
						Indices.rbegin()->second.reserve(WELL_PACKED_VECTOR_SIZE);
						currentIndexPtr = &Indices.rbegin()->second;
						vectorPosition = 0;
					}
				}

				if (1) {
					const char* currentCharPtr = line.c_str();
					while (*currentCharPtr != 0) {
						char code = base2Code(*currentCharPtr);
						(*currentIndexPtr)[vectorPosition++] = (bool)((code << 5) >> 7);
						(*currentIndexPtr)[vectorPosition++] = (bool)((code << 6) >> 7);
						(*currentIndexPtr)[vectorPosition++] = (bool)((code << 7) >> 7);
						++currentCharPtr;
					}
				}
			} // iter file
			{
				int basesInChunk = vectorPosition / 3;
				int endBaseIndex =  Indices.rbegin()->first.StartBP + basesInChunk;
				Indices.rbegin()->first.EndBP = endBaseIndex;
			}
		} // iter chromosomes
		cout << "loaded, " << Indices.size() << "\n";
		cout << "created\n";
	}


	string Search(const TChromInterval& searchInterval) const {
		string output;
		for (int chunkIndex = 0; chunkIndex < Indices.size(); ++chunkIndex) {
			const TChromInterval& interval = Indices[chunkIndex].first;
			if (interval.ChromCode != searchInterval.ChromCode || interval.EndBP <= searchInterval.StartBP ) {
				continue;
			}
			if (interval.StartBP >= searchInterval.EndBP) {
				break;
			}
			//cout << Indices[chunkIndex].second[0] << Indices[chunkIndex].second[1] << Indices[chunkIndex].second[2] << "\n";
			int start = std::max(interval.StartBP, searchInterval.StartBP) - interval.StartBP;
			int end = std::min(interval.EndBP, searchInterval.EndBP) - interval.StartBP;
			const vector<bool>& values = Indices[chunkIndex].second;
			for (int position = start; position < end; ++position) {
				int start = position * 3;
				bool first = values.at(start);
				bool second = values.at(start + 1);
				bool third = values.at(start + 2);
				output += BASE_CODES[(first << 2) + (second << 1) + third];
			}
		}
		return output;
	}


private:
	vector<pair<TChromInterval, vector<bool> > > Indices;
};



TSequenceDB db;


void http::server::request_handler::handle_request(const request& req, reply& rep) {
	{
		bool right_prefix = req.uri.find(SEQUENCE_SERVICE_PREFIX, 0) == 0;
		if (!right_prefix) {
		    rep = reply::stock_reply(reply::bad_request);
		    return;
		}
	}
	vector<TChromInterval> queries = parseQuery(req);
	string responseString;
	for (int queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
		responseString.append(db.Search(queries[queryIndex]));
		responseString.append("\n");
		responseString += "=====\n";
	}
	// Fill out the reply to be sent to the client.
	rep.status = reply::ok;
	rep.content = responseString;
	rep.headers.resize(2);
	rep.headers[0].name = "Content-Length";
	rep.headers[0].value = boost::lexical_cast<std::string>(rep.content.size());
	rep.headers[1].name = "Content-Type";
	rep.headers[1].value = mime_types::extension_to_type("");
}



int main(int argc, char* argv[]) {
	try {
		if (argc != 3) {
			std::cerr << "Usage: http_server <port> <doc_root>\n";
			return 1;
		}
		// Initialise the server.
		http::server::server server("127.0.0.1", argv[1], argv[2]);
		server.run();
	} catch (std::exception& e) {
		std::cerr << "exception: " << e.what() << "\n";
	}

	return 0;
}
