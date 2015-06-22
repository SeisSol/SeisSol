/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
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
 */

#include "vertex_hashmap.h"
#include <iostream>
extern "C"{
	void hash_map_create_instance_(long long& ref){
		ref=(long long) new VertexHashmap();
	}

	void hash_map_destroy_instance_(long long& ref){
		delete (VertexHashmap*) ref;
	}

	void hash_map_get_vertex_id_(long long& ref,double &x,double &y,double &z,int& id){
		((VertexHashmap*) ref)->get(x,y,z,id);

	}
	void hash_map_get_vid_(long long& ref,int& vid,int& id){
		((VertexHashmap*) ref)->get_vid(vid,id);
	}
}

VertexHashmap::VertexHashmap():
	_vertexCounter(0){

}

VertexHashmap::~VertexHashmap(){

}

void VertexHashmap::get(const double x,const double y,const double z,int& id){
	std::vector<double> key(3,0.0);

	key[0]=x;
	key[1]=y;
	key[2]=z;

	std::map<std::vector<double>,int>::iterator it =_map.find(key);
	if(it!=_map.end()){
		_vids.push_back((*it).second);
		id=-1;
	}else{
		_map[key]=_vertexCounter;
		_vids.push_back(_vertexCounter);
		id=_vertexCounter++;
	}

}

void VertexHashmap::get_vid(const int vid,int& id){
	id=_vids[vid];
}
