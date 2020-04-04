classdef feSolve
    %feSolve Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Definition
        Mesh feMesh
        ProblemParameters
        SolverParameters feSolveParameters
        
        % Solver Status
        step_time
    end
    
    methods
        function obj = feSolve(Mesh, ProblemParameters, SolverParameters)
            %feSolve Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 2
                obj.Mesh = Mesh;
                obj.ProblemParameters = ProblemParameters;
                obj.SolverParameters = feSolveParameters();
            elseif nargin == 3
                obj.Mesh = Mesh;
                obj.ProblemParameters = ProblemParameters;
                obj.SolverParameters = SolverParameters;
            end
        end
    end
    
    methods
        % Solvers
        function sol = solve(obj)
            %SOLVE Summary of this method goes here
            %   Detailed explanation goes here
            num_elems = length(obj.Mesh.Elements);
            
            curr_step = 0;
            max_step = obj.SolverParameters.MaxSteps;
            curr_time = 0;
            max_time = 1;
            time = linspace(curr_time, max_time, max_step);
            
            while curr_step < max_step && curr_time < max_time
                curr_step = curr_step+1;
                obj.step_time = time(curr_step);
                
                % Initialize Initial Guess
                for e = 1:num_elems
                    num_elem_nodes = length(obj.Mesh.Elements(e).NodeConnectivity);
                    for n = 1:num_elem_nodes
                        num_node_dof = length(obj.Mesh.Elements(e).Reference.DirichletConditions{n});
                        for d = 1:num_node_dof
                            dof_value = obj.Mesh.Elements(e).Reference.DirichletConditions{n}(d);
                            if isnan(dof_value) == false
                                % Apply scaled boundary condition as solution
                                scaled_dof_value = obj.step_time * dof_value;
                                obj.Mesh.Elements(e).Reference.Solution{n,1}(d,1) = scaled_dof_value;
                            elseif curr_step == 1
                                % Initialize DOF to zero
                                obj.Mesh.Elements(e).Reference.Solution{n,1}(d,1) = 0;
                            end
                        end
                    end
                end
                
                % Perform Nonlinear Newton-Raphson Solver
                sol = newton(obj);
            end
        end
        
        function obj = newton(obj)
            num_elems = length(obj.Mesh.Elements);

            tol = obj.SolverParameters.RelNonlinearTolerance;
            max_iter = obj.SolverParameters.MaxNonlinearIterations;
            curr_iter = 0;
            while curr_iter < max_iter
                curr_iter = curr_iter + 1;
                
                % Evaluate material model
                for e = 1:num_elems
                    obj.Mesh.Elements(e).Reference.MaterialConstitutiveMatrix = feElement.compute_material_consitutive_matrix(1,0);
                end
                
                % Evaluate local stiffness matrices
                for e = 1:num_elems
                obj.Mesh.Elements(e).Reference.StiffnessMatrix = obj.Mesh.Elements(e).compute_local_stiffness_matrix();
                end
                
                % Evaluate local external force vector
                for e = 1:num_elems
                    obj.Mesh.Elements(e).Reference.ExternalForceVector = obj.Mesh.Elements(e).compute_local_externalforce_vector();
                end
                
                % Evaluate local internal force vector
                for e = 1:num_elems
                    d_coeff = obj.Mesh.Elements(e).Reference.Solution;
                    obj.Mesh.Elements(e).Reference.InternalForceVector = obj.Mesh.Elements(e).compute_local_internalforce_vector(d_coeff);
                end
                
                % Assemble global systems
                K = obj.Mesh.assemble_global_stiffness_matrix();
                Fext = obj.Mesh.assemble_global_externalforce_vector();
                Fint = obj.Mesh.assemble_global_internalforce_vector();
                Res = Fext - Fint;
                              
                % Apply Dirichlet BCs to global system
                % Grab all constrained dofs from the elements
                const_dof_list = [];
                const_dof_value = [];
                for e = 1:num_elems
                    num_elem_nodes = length(obj.Mesh.Elements(e).NodeConnectivity);
                    for n = 1:num_elem_nodes
                        num_node_dofs = size(obj.Mesh.Elements(e).DOFConnectivity,1);
                        for d = 1:num_node_dofs
                            dof_value = obj.Mesh.Elements(e).Reference.DirichletConditions{n}(d);
                            if isnan(dof_value) == false
                                const_dof_list = [const_dof_list; obj.Mesh.Elements(e).DOFConnectivity(d,n)];
                                const_dof_value = [const_dof_value; dof_value];
                            end
                        end
                    end
                end
                
                % Condense to just the unique dof entries
                [~, sort_index, ~] = unique(const_dof_list);
                const_dof_list = const_dof_list(sort_index);
                const_dof_value = const_dof_value(sort_index);
                
                % Subtract known equations from RHS
                Res(const_dof_list) = Res(const_dof_list) - (K(const_dof_list,const_dof_list) * const_dof_value);
                
                % Remove known equations from system of equations
                K(const_dof_list,:) = [];
                K(:,const_dof_list) = [];
                Res(const_dof_list) = [];
                
                % Test for Convergence
                if norm(Res) < tol
                    disp("Converged!")
                    return
                end
                
                % Solve Linear Solution
                du = K\Res;
                
                % Reassemble Solution
                dof_list = (1:max(obj.Mesh.DOFConnectivity(:)))';
                dof_list(const_dof_list) = [];
                for e = 1:num_elems
                    num_elem_nodes = length(obj.Mesh.Elements(e).NodeConnectivity);
                    for n = 1:num_elem_nodes
                        num_node_dofs = size(obj.Mesh.Elements(e).DOFConnectivity,1);
                        for node_dof = 1:num_node_dofs
                            for d = 1:length(dof_list)
                                if dof_list(d) == obj.Mesh.Elements(e).DOFConnectivity(node_dof,n)
                                    obj.Mesh.Elements(e).Reference.Solution{n}(node_dof) = obj.Mesh.Elements(e).Reference.Solution{n}(node_dof) + du(node_dof);
                                end
                            end
                        end
                    end
                end                
            end
            
        end
    end
end

