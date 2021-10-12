Fault tagging
=============

| In SeisSol, boundary conditions are tagged as:
| 0: regular
| 1: free surface
| 2: free surface + gravity (water surface)
| 3 or n>64: dynamic rupture
| 5: absorbing
| 6: periodic

Dynamic rupture can therefore be tagged using a range of possible tags.
This allows initializing fault parameters segment-wise
easily. For example, if we have 2 segments, and we want them to have
different dynamic friction, we can tag them with 3 and 65 and then use:

.. code-block:: yaml

   [mu_d]: !Any
     components:
       - !GroupFilter
         groups: 3
         components: !ConstantMap
           map:
             mu_d:    0.3
       - !GroupFilter
         groups: 65
         components: !ConstantMap
           map:
             mu_d:    0.4

Currently, the only way to tag fault faces other tags than 3 with SimModeler is to use the `--xml` option of pumgen. 
For example, to tag face 2 as 3 and face 8 and 9 as 65, we would
use:

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition>
   <boundaryCondition tag="65">8,9</boundaryCondition>

Then pumgen is run using the xml option:

::

   pumgen -s simmodsuite -l SimModelerLib.lic --xml MeshandAnalysisAttributes.xml prefix.smd output_prefix


Using more than 189 dynamic rupture tags
----------------------------------------

Currently, SeisSol cannot handle more than 255 fault tags, that is 189 dynamic rupture tags. To overcome this limitation, it is necessary to patch PUMGen, SeisSol and the PUML submodule of SeisSol. This can be done with:

- Download all 3 files containing necessary changes.
- patch_PUMGen.diff (https://github.com/palgunadi1993/SeisSol-PUMGen-PUML/blob/master/patch_PUMGen.diff)

.. code:: diff

    diff --git a/src/pumgen.cpp b/src/pumgen.cpp
    index 5a19c60..794f5bb 100644
    --- a/src/pumgen.cpp
    +++ b/src/pumgen.cpp
    @@ -332,7 +332,7 @@ int main(int argc, char* argv[]) {
         checkH5Err(h5space);
     
         hid_t h5group =
    -        H5Dcreate(h5file, "/group", H5T_STD_I32LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    +        H5Dcreate(h5file, "/group", H5T_STD_I64LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         checkH5Err(h5group);
     
         start[0] = offsets[0];
    @@ -343,18 +343,21 @@ int main(int argc, char* argv[]) {
         h5memspace = H5Screate_simple(1, sizes, 0L);
         checkH5Err(h5memspace);
     
    -    int* group = new int[localSize[0]];
    +    int64_t* group = new int64_t[localSize[0]];
         it = mesh->begin(3);
         index = 0;
         while (apf::MeshEntity* element = mesh->iterate(it)) {
           assert(mesh->hasTag(element, groupTag));
    +      
    +      int myGroup; 
    +      mesh->getIntTag(element, groupTag, &myGroup);
    +      group[index] = static_cast<int64_t>(myGroup);
     
    -      mesh->getIntTag(element, groupTag, &group[index]);
           index++;
         }
         mesh->end(it);
     
    -    checkH5Err(H5Dwrite(h5group, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, group));
    +    checkH5Err(H5Dwrite(h5group, H5T_NATIVE_INT64, h5memspace, h5space, h5dxlist, group));
     
         checkH5Err(H5Sclose(h5space));
         checkH5Err(H5Sclose(h5memspace));
    @@ -370,15 +373,15 @@ int main(int argc, char* argv[]) {
       apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
       assert(boundaryTag);
     
    -  int* boundary = new int[localSize[0]];
    -  memset(boundary, 0, localSize[0] * sizeof(int));
    +  int64_t* boundary = new int64_t[localSize[0]];
    +  memset(boundary, 0, localSize[0] * sizeof(int64_t));
     
       sizes[0] = globalSize[0];
       h5space = H5Screate_simple(1, sizes, 0L);
       checkH5Err(h5space);
     
    -  hid_t h5boundary =
    -      H5Dcreate(h5file, "/boundary", H5T_STD_I32LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    +  hid_t h5boundary = 
    +      H5Dcreate(h5file, "/boundary", H5T_STD_I64LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       checkH5Err(h5boundary);
     
       start[0] = offsets[0];
    @@ -397,13 +400,15 @@ int main(int argc, char* argv[]) {
     
         for (unsigned int i = 0; i < 4; i++) {
           if (mesh->hasTag(faces[i], boundaryTag)) {
    -        int b;
    -        mesh->getIntTag(faces[i], boundaryTag, &b);
    +        int ba;
    +        int64_t b;
    +        mesh->getIntTag(faces[i], boundaryTag, &ba);
    +        b = static_cast<int64_t>(ba);
     
    -        if (b <= 0 || b > std::numeric_limits<unsigned char>::max())
    +        if (b <= 0 || b > std::numeric_limits<unsigned short>::max())
               logError() << "Cannot handle boundary condition" << b;
     
    -        boundary[index] += b << (i * 8);
    +        boundary[index] += b << (i * 16);
           }
         }
     
    @@ -411,7 +416,7 @@ int main(int argc, char* argv[]) {
       }
       mesh->end(it);
     
    -  checkH5Err(H5Dwrite(h5boundary, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, boundary));
    +  checkH5Err(H5Dwrite(h5boundary, H5T_NATIVE_INT64, h5memspace, h5space, h5dxlist, boundary));
     
       checkH5Err(H5Sclose(h5space));
       checkH5Err(H5Sclose(h5memspace));


- patch_SeisSol.diff (https://github.com/palgunadi1993/SeisSol-PUMGen-PUML/blob/master/patch_SeisSol.diff)

.. code:: diff

    diff --git a/src/Geometry/PUMLReader.cpp b/src/Geometry/PUMLReader.cpp
    index 23e5379a..ba0115f7 100644
    --- a/src/Geometry/PUMLReader.cpp
    +++ b/src/Geometry/PUMLReader.cpp
    @@ -313,8 +313,8 @@ void seissol::PUMLReader::getMesh(const PUML::TETPUML &puml)
     	const std::vector<PUML::TETPUML::face_t> &faces = puml.faces();
     	const std::vector<PUML::TETPUML::vertex_t> &vertices = puml.vertices();
     
    -	const int* material = puml.cellData(0);
    -	const int* boundaryCond = puml.cellData(1);
    +	const int64_t* material = puml.cellData(0);
    +	const int64_t* boundaryCond = puml.cellData(1);
     
     	std::unordered_map<int, std::vector<unsigned int> > neighborInfo; // List of shared local face ids
     
    @@ -369,8 +369,8 @@ void seissol::PUMLReader::getMesh(const PUML::TETPUML &puml)
     				m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
     			}
     
    -			int bcCurrentFace = (boundaryCond[i] >> (j*8)) & 0xFF;
    -			int faultTag = bcCurrentFace;
    +			int64_t bcCurrentFace = (boundaryCond[i] >> (j*16)) & 0xFFFF;
    +			int64_t faultTag = bcCurrentFace;
     			if (bcCurrentFace > 64) {
     				bcCurrentFace = 3;
     			}
    @@ -379,7 +379,7 @@ void seissol::PUMLReader::getMesh(const PUML::TETPUML &puml)
     			m_elements[i].mpiIndices[FACE_PUML2SEISSOL[j]] = 0;
     		}
     
    -		m_elements[i].material = material[i];
    +		m_elements[i].material = static_cast<int>(material[i]);
     	}
     
     	// Exchange ghost layer information and generate neighbor list
    diff --git a/src/Initializer/ParameterDB.cpp b/src/Initializer/ParameterDB.cpp
    index 18529caa..326b95d1 100644
    --- a/src/Initializer/ParameterDB.cpp
    +++ b/src/Initializer/ParameterDB.cpp
    @@ -85,7 +85,7 @@ easi::Query seissol::initializers::ElementBarycentreGeneratorPUML::generate() co
       std::vector<PUML::TETPUML::cell_t> const& cells = m_mesh.cells();
       std::vector<PUML::TETPUML::vertex_t> const& vertices = m_mesh.vertices();
     
    -  int const* material = m_mesh.cellData(0);
    +  int64_t const* material = m_mesh.cellData(0);
       
       easi::Query query(cells.size(), 3);
       for (unsigned cell = 0; cell < cells.size(); ++cell) {
    @@ -105,7 +105,7 @@ easi::Query seissol::initializers::ElementBarycentreGeneratorPUML::generate() co
           query.x(cell,dim) *= 0.25;
         }
         // Group
    -    query.group(cell) = material[cell];
    +    query.group(cell) = static_cast<int>(material[cell]);
       }
       return query;
     }
    diff --git a/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp b/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp
    index ddc9b856..cfcd6b67 100644
    --- a/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp
    +++ b/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp
    @@ -152,8 +152,8 @@ int LtsWeights::getCluster(double timestep, double globalMinTimestep, unsigned r
       return cluster;
     }
     
    -int LtsWeights::getBoundaryCondition(int const *boundaryCond, unsigned cell, unsigned face) {
    -  int bcCurrentFace = ((boundaryCond[cell] >> (face * 8)) & 0xFF);
    +int LtsWeights::getBoundaryCondition(int64_t const *boundaryCond, unsigned cell, unsigned face) {
    +  int bcCurrentFace = ((boundaryCond[cell] >> (face * 16)) & 0xFFFF);
       if (bcCurrentFace > 64) {
         bcCurrentFace = 3;
       }
    @@ -227,7 +227,7 @@ std::vector<int> LtsWeights::computeCostsPerTimestep() {
       const auto &cells = m_mesh->cells();
     
       std::vector<int> cellCosts(cells.size());
    -  int const *boundaryCond = m_mesh->cellData(1);
    +  int64_t const *boundaryCond = m_mesh->cellData(1);
       for (unsigned cell = 0; cell < cells.size(); ++cell) {
         int dynamicRupture = 0;
         int freeSurfaceWithGravity = 0;
    @@ -269,7 +269,7 @@ int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
     
       std::vector<PUML::TETPUML::cell_t> const &cells = m_mesh->cells();
       std::vector<PUML::TETPUML::face_t> const &faces = m_mesh->faces();
    -  int const *boundaryCond = m_mesh->cellData(1);
    +  int64_t const *boundaryCond = m_mesh->cellData(1);
     
     #ifdef USE_MPI
       std::unordered_map<int, std::vector<int>> rankToSharedFaces;
    @@ -377,4 +377,4 @@ int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
     
       return numberOfReductions;
     }
    -} // namespace seissol::initializers::time_stepping
    \ No newline at end of file
    +} // namespace seissol::initializers::time_stepping
    diff --git a/src/Initializer/time_stepping/LtsWeights/LtsWeights.h b/src/Initializer/time_stepping/LtsWeights/LtsWeights.h
    index a3c65b5d..c3ebf0ab 100644
    --- a/src/Initializer/time_stepping/LtsWeights/LtsWeights.h
    +++ b/src/Initializer/time_stepping/LtsWeights/LtsWeights.h
    @@ -86,7 +86,7 @@ protected:
       GlobalTimeStepDetails collectGlobalTimeStepDetails();
       void computeMaxTimesteps(std::vector<double> const &pWaveVel, std::vector<double> &timeSteps);
       int getCluster(double timestep, double globalMinTimestep, unsigned rate);
    -  int getBoundaryCondition(int const *boundaryCond, unsigned cell, unsigned face);
    +  int getBoundaryCondition(int64_t const *boundaryCond, unsigned cell, unsigned face);
       std::vector<int> computeClusterIds();
       int enforceMaximumDifference();
       int enforceMaximumDifferenceLocal(int maxDifference = 1);


- patch_PUML.diff (https://github.com/palgunadi1993/SeisSol-PUMGen-PUML/blob/master/patch_PUML.diff)

.. code:: diff

    diff --git a/PUML.h b/PUML.h
    index b0735be..802b3d8 100644
    --- a/PUML.h
    +++ b/PUML.h
    @@ -117,7 +117,7 @@ private:
     	internal::VertexElementMap<2> m_v2e;
     
     	/** User cell data */
    -	std::vector<int*> m_cellData;
    +	std::vector<int64_t*> m_cellData;
     
     	/** User vertex data */
     	std::vector<int*> m_vertexData;
    @@ -138,7 +138,7 @@ public:
     		delete [] m_originalCells;
     		delete [] m_originalVertices;
     
    -		for (std::vector<int*>::const_iterator it = m_cellData.begin();
    +		for (std::vector<int64_t*>::const_iterator it = m_cellData.begin();
     				it != m_cellData.end(); ++it) {
     			delete [] *it;
     		}
    @@ -338,8 +338,8 @@ public:
     		checkH5Err(H5Pset_dxpl_mpio(h5alist, H5FD_MPIO_COLLECTIVE));
     #endif // USE_MPI
     
    -		int* data = new int[localSize];
    -		checkH5Err(H5Dread(h5dataset, H5T_NATIVE_INT, h5memspace, h5space, h5alist, data));
    +		int64_t* data = new int64_t[localSize];
    +		checkH5Err(H5Dread(h5dataset, H5T_NATIVE_INT64, h5memspace, h5space, h5alist, data));
     
     		// Close data
     		checkH5Err(H5Sclose(h5space));
    @@ -356,7 +356,8 @@ public:
     			m_cellData.push_back(data);
     			break;
     		case VERTEX:
    -			m_originalVertexData.push_back(data);
    +            int* aa = (int*)data;
    +			m_originalVertexData.push_back(aa);
     			break;
     		}
     	}
    @@ -398,9 +399,9 @@ public:
     		m_originalCells = newCells;
     
     		// Sort other data
    -		for (std::vector<int*>::iterator it = m_cellData.begin();
    +		for (std::vector<int64_t*>::iterator it = m_cellData.begin();
     				it != m_cellData.end(); ++it) {
    -			int* newData = new int[m_originalSize[0]];
    +			int64_t* newData = new int64_t[m_originalSize[0]];
     			for (unsigned int i = 0; i < m_originalSize[0]; i++) {
     				newData[i] = (*it)[indices[i]];
     			}
    @@ -455,11 +456,11 @@ public:
     		MPI_Type_free(&cellType);
     
     		// Exchange cell data
    -		for (std::vector<int*>::iterator it = m_cellData.begin();
    +		for (std::vector<int64_t*>::iterator it = m_cellData.begin();
     				it != m_cellData.end(); ++it) {
    -			int* newData = new int[m_originalSize[0]];
    -			MPI_Alltoallv(*it, sendCount, sDispls, MPI_INT,
    -				newData, recvCount, rDispls, MPI_INT,
    +			int64_t* newData = new int64_t[m_originalSize[0]];
    +			MPI_Alltoallv(*it, sendCount, sDispls, MPI_INT64_T,
    +				newData, recvCount, rDispls, MPI_INT64_T,
     				m_comm);
     
     			delete [] *it;
    @@ -907,7 +908,7 @@ public:
     	/**
     	 * @return User cell data
     	 */
    -	const int* cellData(unsigned int index) const
    +	const int64_t* cellData(unsigned int index) const
     	{
     		return m_cellData[index];
     	}


Then apply the patches:

.. code-block::
  cd PUMGen
  git apply patch_PUMGen.diff


.. code-block::
  cd SeisSol
  git apply patch_SeisSol.diff


and finally:

.. code-block::
  cd SeisSol/submodules/PUML/
  git apply patch_PUML.diff


Meshes with more than 255 tags can be created using pumgen -xml option, e.g. :

.. code-block:: xml
   <boundaryCondition tag="3">13245</boundaryCondition>
   .
   .
   .
   <boundaryCondition tag="900">12345,14325</boundaryCondition>


