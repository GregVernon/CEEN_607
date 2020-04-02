clear
M = importMesh("mesh/test_1x1.mat");

num_elems = length(M.Elements);
e = 1;

d_coeff = {[0; 0];[0.1; 0];[0; 0];[0.1; 0]};


M.Elements(e).Reference.MaterialConstitutiveMatrix = feElement.compute_material_consitutive_matrix(1,0);
M.Elements(e).Reference.DirichletConditions = {[0 0]; [0.1 missing]; [0 0]; [0.1 missing]};
M.Elements(e).Reference.BodyForce = [0; 0];
M.Elements(e).Reference.NodeForce = {[0; 0]; [0; 0]; [0; 0]; [0; 0]};
M.Elements(e).Reference.SurfacePressure = {[0]; [0]; [0]; [0]};

M.Elements(e).Reference.StiffnessMatrix = M.Elements(e).compute_local_stiffness_matrix();
M.Elements(e).Reference.InternalForceVector = M.Elements(e).compute_local_internalforce_vector(d_coeff);
M.Elements(e).Reference.ExternalForceVector = M.Elements(e).compute_local_externalforce_vector();
