#include <iostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <fstream>
#include <set>
#include <sstream>
#include "../server/server.hpp"
#include <algorithm>



/*
namespace http {
namespace server {
void request_handler::handle_request(const request& req, reply& rep)
{
  // Decode url to path.
  std::string request_path;
  if (!url_decode(req.uri, request_path))
  {
    rep = reply::stock_reply(reply::bad_request);
    return;
  }

  // Request path must be absolute and not contain "..".
  if (request_path.empty() || request_path[0] != '/'
      || request_path.find("..") != std::string::npos)
  {
    rep = reply::stock_reply(reply::bad_request);
    return;
  }

  // If path ends in slash (i.e. is a directory) then add "index.html".
  if (request_path[request_path.size() - 1] == '/')
  {
    request_path += "index.html";
  }

  // Determine the file extension.
  std::size_t last_slash_pos = request_path.find_last_of("/");
  std::size_t last_dot_pos = request_path.find_last_of(".");
  std::string extension;
  if (last_dot_pos != std::string::npos && last_dot_pos > last_slash_pos)
  {
    extension = request_path.substr(last_dot_pos + 1);
  }

  // Open the file to send back.
  std::string full_path = doc_root_ + request_path;
  std::ifstream is(full_path.c_str(), std::ios::in | std::ios::binary);
  if (!is)
  {
    rep = reply::stock_reply(reply::not_found);
    return;
  }

  // Fill out the reply to be sent to the client.
  rep.status = reply::ok;
  char buf[512];
  while (is.read(buf, sizeof(buf)).gcount() > 0)
    rep.content.append(buf, is.gcount());
  rep.headers.resize(2);
  rep.headers[0].name = "Content-Length";
  rep.headers[0].value = boost::lexical_cast<std::string>(rep.content.size());
  rep.headers[1].name = "Content-Type";
  rep.headers[1].value = mime_types::extension_to_type(extension);
}

}
}
*/

using std::pair;
using std::string;
using std::vector;
using std::cout;

const string SNP_SERVICE_PREFIX = "/js_snp?";
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

struct TQueryInterval {
	int ChromCode;
	int StartBP, EndBP;
	TQueryInterval(const string& chromosome, const int start, const int end) :
		ChromCode(chrom2chromCode(chromosome)),
		StartBP(start),
		EndBP(end) {
	}
	TQueryInterval() {
	}
};




vector<TQueryInterval> parseQuery(const http::server::request& req) {
	vector<TQueryInterval> queries;
	{//parse
		string query = req.uri.substr(SNP_SERVICE_PREFIX.size(), req.uri.size() - SNP_SERVICE_PREFIX.size());
		vector<string> chunks = splitString(query, QUERY_DELIM);
		for (vector<string>::const_iterator chunkIt = chunks.begin(); chunkIt != chunks.end(); ++chunkIt) {
			vector<string> queryFields = splitString(*chunkIt, FIELD_DELIM);
			try {
				string queryChromosome = queryFields[0];
				int startBP = boost::lexical_cast<int>(queryFields[1]);
				int endBP = boost::lexical_cast<int>(queryFields[2]);
				queries.push_back(TQueryInterval(queryChromosome, startBP, endBP));
			} catch (std::exception& e) {
				std::cerr << "exception: " << e.what() << "\n";
			}
		}
	}
	return queries;
}


const char* DATA_FILE = "../data/snp/all_rs.txt";
const string DATA_FILE_FIELD_DELIM = " ";




struct TSNP {
	int BpIndex;
	int ID;
	string FirstAllele, SecondAllele;
	TSNP(int bpIndex, int id, const string& firstAllele, const string& secondAllele) :
		BpIndex(bpIndex), ID(id), FirstAllele(firstAllele), SecondAllele(secondAllele) {

	}
	bool operator>(const TSNP& other) const {
		return BpIndex > other.BpIndex;
	}
	bool operator<(const TSNP& other) const {
		return BpIndex < other.BpIndex;
	}
	bool operator==(const TSNP& other) const {
		return BpIndex == other.BpIndex;
	}
};


template <class TKey, class TValue>
class TBinarySearchMap {
public:
	typedef typename vector<pair<TKey, TValue> >::const_iterator TResponseIterator;

	TBinarySearchMap() {
	}
	//destroys other!!!
	void operator=(TBinarySearchMap& other) {
		Data.swap(other.Data);
	}

	//data vector will be modified
	TBinarySearchMap(vector<pair<TKey, TValue> >* dataPtr) {
		Data.swap(*dataPtr);
		std::sort(Data.begin(), Data.end());
	}

	TResponseIterator FindFirst(const TKey& key) const {
		if (Data.empty()) {
			return End();
		}
		int left = 0;
		int right = (int)Data.size() - 1;
		if (Data[right].first < key) {
			return End();
		}
		if (Data[left].first >= key) {
			return Data.begin();
		}
		while (right - left > 1) {
			int middle = (left + right) >> 1;
			if (Data[middle].first < key) {
				left = middle;
			} else {
				right = middle;
			}
		}
		return Data.begin() + right;
	}

	TResponseIterator End() const {
		return Data.end();
	}

private:
	vector<pair<TKey, TValue> > Data;

};


const int WELL_PACKED_VECTOR_SIZE = 16777216;
bool DEBUG = true;

class TSNP_DB {
public:
	TSNP_DB() {
		cout << "started\n";
		vector<pair<TKey, TSNP> > snps;
		{
			std::ifstream file(DATA_FILE);
			if (!file.is_open()) {
				return;
			}

			string line;
			char buffer[100];
			while (getline(file, line)) {
				if (line.c_str()[0] == '#') {
					continue;
				}
				if (line.find("rs") == line.npos) {
					continue;
				}
				if (snps.size() % 1000000 == 0) {
					cout << "..processed: " << snps.size() << "\n";
				}
				//I want it to be fast
				int first_delim = line.find(DATA_FILE_FIELD_DELIM);
				int second_delim = line.find(DATA_FILE_FIELD_DELIM, first_delim + 1);
				int third_delim = line.find(DATA_FILE_FIELD_DELIM, second_delim + 1);
				int fourth_delim = line.find(DATA_FILE_FIELD_DELIM, third_delim + 1);
				if (first_delim == string::npos || second_delim == string::npos ||
						third_delim == string::npos || fourth_delim == string::npos) {
					cout << "FUCKUP " << line << "\n";
					exit(0);
				}
				int chromCode = chrom2chromCode(line);
				int bpIndex = 0;
				{
					for (int pos = first_delim + 1; pos < second_delim; ++pos) {
						buffer[pos - first_delim - 1] = line.c_str()[pos];
					}
					buffer[second_delim - first_delim - 1] = 0;
					bpIndex = atoi((const char*)buffer);
				}
				int rsId = 0;
				{
					for (int pos = second_delim + 3; pos < third_delim; ++pos) {
						buffer[pos - second_delim - 3] = line.c_str()[pos];
					}
					buffer[third_delim - second_delim - 3] = 0;
					rsId = atoi((const char*)buffer);
				}
				string firstAllele = line.substr(third_delim + 1, fourth_delim - third_delim - 1);
				string secondAllele = line.substr(fourth_delim + 1, (int)line.size() - fourth_delim - 1);

				//if (chromCode == 1 && bpIndex < 1000000 && bpIndex > 998000) {
				//	cout << "data: " << bpIndex << ", " << rsId << ", " << firstAllele << "\n";
				//}

				if (snps.size() >= WELL_PACKED_VECTOR_SIZE) {
					TBinarySearchMap<TKey, TSNP> index(&snps);
					Indices.push_back(TBinarySearchMap<TKey, TSNP>());
					*Indices.rbegin() = index;
					snps.clear();
					if (DEBUG) break;
				}
				snps.push_back(pair<TKey, TSNP>(TKey(chromCode, bpIndex), TSNP(bpIndex, rsId, firstAllele, secondAllele)));
			}
			Indices.push_back( TBinarySearchMap<TKey, TSNP>(&snps));
		}
		cout << "loaded, " << snps.size() << "\n";
		cout << "created\n";
	}

	vector<TSNP> Search(const TQueryInterval& query) {
		TKey startKey = TKey(query.ChromCode, query.StartBP);
		TKey endKey = TKey(query.ChromCode, query.EndBP);
		vector<TSNP> response;
		for (int index = 0; index < Indices.size(); ++index) {
			TBinarySearchMap<TKey, TSNP>::TResponseIterator iter = Indices[index].FindFirst(startKey);
			while (iter != Indices[index].End() && iter->first < endKey) {
				response.push_back(iter->second);
				//cout << iter->second.BpIndex << " " << iter->second.ID << " " << iter->second.SecondAllele << "\n";
				++iter;
			}
		}
		return response;
	}

private:
	typedef pair<char, int> TKey;
	vector<TBinarySearchMap<TKey, TSNP> > Indices;
};



TSNP_DB db;

void http::server::request_handler::handle_request(const request& req, reply& rep) {
	{
		bool right_prefix = req.uri.find(SNP_SERVICE_PREFIX, 0) == 0;
		if (!right_prefix) {
		    rep = reply::stock_reply(reply::bad_request);
		    return;
		}
	}
	vector<TQueryInterval> queries = parseQuery(req);
	string responseString;
	for (int queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
		vector<TSNP> responses = db.Search(queries[queryIndex]);
		for (int respIndex = 0; respIndex < responses.size(); ++respIndex) {
			responseString.append(boost::lexical_cast<std::string>(responses[respIndex].BpIndex));
			responseString.append(" ");
			responseString.append(boost::lexical_cast<std::string>(responses[respIndex].ID));
			responseString.append(" ");
			responseString.append(responses[respIndex].FirstAllele + " " + responses[respIndex].SecondAllele);
			responseString.append("\n");
		}
		responseString += "\n";
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
