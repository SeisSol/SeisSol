/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Parser for XML-files.
 **/

#ifndef XMLPARSER_HPP_
#define XMLPARSER_HPP_

#include <cassert>
#include <string>
#include <pugixml/src/pugixml.hpp>
#include <utils/logger.h>

namespace seissol {
  class XmlParser;
}

class seissol::XmlParser {
  //private:
    //!XML-document, which holds information about the setup of the matrix kernels.
    pugi::xml_document m_matrices;

  public:
    /**
     * Constructor, which loads the matrix-XML-file.
     *
     * @param i_pathToMatrixXmlFile location of the xml-file, which contains the matrix information.
     **/
    XmlParser( std::string i_pathToMatrixXmlFile ) {
     logDebug() << "loading XML-file: " << i_pathToMatrixXmlFile;

     // load the xml-file
     pugi::xml_parse_result l_result = m_matrices.load_file(i_pathToMatrixXmlFile.c_str());
     
     if( !l_result ) {
       logError() << "XML-load failed: " << l_result.description() << ", File: " << i_pathToMatrixXmlFile;
     }
    }

    /**
     * Read global matrices of the specified type.
     *
     * @param i_type matrix type.
     * @param io_ids ids of the matrices.
     * @param io_names names of the matrices.
     * @param io_numberOfRows number of rows.
     * @param io_numberOfColumns number of columns.
     * @param io_sparsities sparsity handling of the matrices (true, if sparse; false, if dense).
     * @param io_rows row indicices of non-zero entries.
     * @param io_columns column indices of non-zero entries.
     * @param io_values values of non-zero entries.
     **/
    void readGlobalMatrices( const std::string                               &i_type,
                                   std::vector< unsigned int >              &io_ids,
                                   std::vector< std::string  >              &io_names,
                                   std::vector< unsigned int >              &io_numberOfRows,
                                   std::vector< unsigned int >              &io_numberOfColumns,
                                   std::vector< bool         >              &io_sparsities,
                                   std::vector< std::vector<unsigned int> > &io_rows,
                                   std::vector< std::vector<unsigned int> > &io_columns,
                                   std::vector< std::vector<double>       > &io_values ) const {
                                     
      pugi::xml_node l_matrix = m_matrices.child("matrices").child("global").child(i_type.c_str());
      if (l_matrix.empty()) {
        logError() << "The XML-File containing the constant matrices of the dG method lacks matrices of the type" << i_type.c_str();
      }
                                     
      for (; l_matrix; l_matrix = l_matrix.next_sibling(i_type.c_str()) ) {
        // read meta-information
        io_ids.push_back(             l_matrix.attribute("id").as_uint()      );
        io_names.push_back(           l_matrix.attribute("name").as_string()  );
        // TODO: reintroduce sparse-dense swich
        io_sparsities.push_back( false );
        //io_sparsities.push_back(      l_matrix.attribute("sparse").as_bool()  );
        io_numberOfRows.push_back(    l_matrix.attribute("rows").as_uint()    );
        io_numberOfColumns.push_back( l_matrix.attribute("columns").as_uint() );

        // read entries
        io_rows.push_back(    std::vector<unsigned int>() );
        io_columns.push_back( std::vector<unsigned int>() );
        io_values.push_back(  std::vector<double>()       );
        for( pugi::xml_node l_entry = l_matrix.child("entry");
             l_entry;
             l_entry = l_entry.next_sibling("entry") ) {
          io_rows.back().push_back(    l_entry.attribute("row").as_uint()      );
          io_columns.back().push_back( l_entry.attribute("column").as_uint()   );
          io_values.back().push_back(  l_entry.attribute("value").as_double()  );
        }
      }
    }

    /**
     * Read the structure of the local matrices of the specified type.
     *
     * @param i_type matrix type.
     * @param io_ids ids of the matrices.
     * @param io_numberOfRows number of rows.
     * @param io_numberOfColumns number of columns.
     * @param io_names names of the matrices.
     * @param io_sparsities sparsity handling of the matrices (true, if sparse; false, if dense).
     * @param io_rows row indicices of non-zero entries.
     * @param io_columns column indices of non-zero entries.
     **/
   void readLocalMatrixStructure( const std::string           &i_type,
                                  std::vector< unsigned int > &io_ids,
                                  std::vector< unsigned int > &io_numberOfRows,
                                  std::vector< unsigned int > &io_numberOfColumns,
                                  std::vector< bool         > &io_sparsities,
                                  std::vector< std::vector<unsigned int> > &io_rows,
                                  std::vector< std::vector<unsigned int> > &io_columns ) const {
      pugi::xml_node l_matrix = m_matrices.child("matrices").child("local").child(i_type.c_str());
      
      // read meta-information
      io_ids.push_back(             l_matrix.attribute("id").as_uint()      );
      // TODO: reintroduce sparse-dense swich
      io_sparsities.push_back( false );
      //io_sparsities.push_back(      l_matrix.attribute("sparse").as_bool()  );
      io_numberOfRows.push_back(    l_matrix.attribute("rows").as_uint()    );
      io_numberOfColumns.push_back( l_matrix.attribute("columns").as_uint() );
        
      // read position of entries
      io_rows.push_back(    std::vector<unsigned int>() );
      io_columns.push_back( std::vector<unsigned int>() );
      for( pugi::xml_node l_entry = l_matrix.child("entry");
           l_entry;
           l_entry = l_entry.next_sibling("entry") ) {
        io_rows.back().push_back(    l_entry.attribute("row").as_uint()      );
        io_columns.back().push_back( l_entry.attribute("column").as_uint()   );
      }
   }
};

#endif
