//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-08
//=======================

#ifndef READ_TXT_TOOLS_H
#define READ_TXT_TOOLS_H

#include <vector>
#include <map>
#include <set>
#include <string>

namespace Read_Txt_Tools
{
	// {
	//  	INPUT_PARAMETERS:       [1.0],
	//  	ecutwfc:        		[200, Ry]
	// }
	std::map<std::string, std::vector<std::string>> read_file_to_map(
		const std::string &file_name,
		const std::set<std::string> &comment_seps);


	// [
	//  	[INPUT_PARAMETERS, 1.0],
	//  	[ecutwfc, 200, Ry]
	// ]
	std::vector<std::vector<std::string>> read_file_to_vector(
		const std::string &file_name,
		const std::set<std::string> &comment_seps);


	// [ecutwfc, 200, Ry]
	std::vector<std::string> split_whitespace(const std::string &str);


	// str:
	//		ecutwfc 200 Ry # XXX
	// return:
	//		ecutwfc 200 Ry
	std::string ignore_comment(const std::string &str, const std::set<std::string> &comment_seps);


	// return_prefix + content_new = content_old
	std::vector<std::vector<std::string>> cut_paragraph(
		std::vector<std::vector<std::string>> &content,
		const std::set<std::string> &labels);


	// data_v:
	//	[
	//		[INPUT_PARAMETERS, 1.0],
	//		[ecutwfc, 200, Ry]
	//	]
	// return:
	//	[INPUT_PARAMETERS, 1.0, ecutwfc, 200, Ry]
	template<typename T>
	std::vector<T> chain(const std::vector<std::vector<T>> &data_v)
	{
		std::vector<T> data_chain;
		for(const std::vector<T> & data_i : data_v)
			data_chain.insert(data_chain.end(), data_i.begin(), data_i.end());
		return data_chain;
	}	
}

#endif