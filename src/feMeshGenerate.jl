include("feDatastructure.jl")
include("feMeshImport.jl")

function generateMesh_2DTensorProduct(xyBounds,NX,NY,degree,outname)
    xmin = xyBounds[1]
    xmax = xyBounds[2]
    ymin = xyBounds[3]
    ymax = xyBounds[4]
    num_elem = NX * NY
    num_nodes = (NX+1) * (NY+1)

    ## Generate mesh 
    proc = run(`py cubitMesh.py --xmin $xmin --xmax $xmax --ymin $ymin --ymax $ymax --nx $NX --ny $NY --outname $outname`)

    ## Import generated mesh
    MeshData = importGenesis("2x2_mesh.g")
    return MeshData
end