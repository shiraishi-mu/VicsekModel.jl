var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = VicsekModel","category":"page"},{"location":"#VicsekModel","page":"Home","title":"VicsekModel","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VicsekModel.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VicsekModel]","category":"page"},{"location":"#VicsekModel.VicsekModelParameters","page":"Home","title":"VicsekModel.VicsekModelParameters","text":"Fields\n\nL: Length of the square box\nρ: Density of particles\nN: Number of particles\nr₀: Interaction radius\nΔt: Time step\nfactor: Factor to determine the initial velocity\nv₀: Initial velocity\nη: Noise strength\nalg: Algorithm to update the position and velocity\n\n\n\n\n\n","category":"type"},{"location":"#VicsekModel.VicsekModelVariables","page":"Home","title":"VicsekModel.VicsekModelVariables","text":"Initialize the variables of the Vicsek model\n\npos: Random position of particles\nvel: Random velocity of particles\nθ: Random angle of particles\nS: Complex number to calculate the average angle of neighbors\n\n\n\n\n\n","category":"type"},{"location":"#VicsekModel.VicsekModelVariables-2","page":"Home","title":"VicsekModel.VicsekModelVariables","text":"Variables of the Vicsek model\n\npos: Position of particles\nvel: Velocity of particles\nθ: Angle of particles\nS: Complex number to calculate the average angle of neighbors\n\n\n\n\n\n","category":"type"},{"location":"#Base.show-Tuple{IO, VicsekModel.VicsekModelParameters}","page":"Home","title":"Base.show","text":"Show the parameters of the Vicsek model\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Tuple{IO, VicsekModel.VicsekModelVariables}","page":"Home","title":"Base.show","text":"Show the variables of the Vicsek model\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.calculate_weights!-Union{Tuple{Tf}, Tuple{NearestNeighbors.NNTree, Matrix{Tf}, Vector{Tf}, Vector{Complex}, Tf, Base.OneTo}} where Tf","page":"Home","title":"VicsekModel.calculate_weights!","text":"Calculate the average angle of neighbors\n\nArguments\n\ntree: KDTree to find the neighbors\npos: Position of particles\nθ: Angle of particles\nS: Complex number to calculate the average angle of neighbors\nr₀: Interaction radius\nix: Index of particles\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.periodic_boundary!-Union{Tuple{Tf}, Tuple{Matrix{Tf}, Tf, Base.OneTo}} where Tf","page":"Home","title":"VicsekModel.periodic_boundary!","text":"Periodic boundary condition\n\nArguments\n\npos: Position of particles\nL: Length of the square box\nix: Index of particles\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.update!-Tuple{Any, Any}","page":"Home","title":"VicsekModel.update!","text":"Holy trait of update function\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.update!-Union{Tuple{Ti}, Tuple{Tf}, Tuple{VicsekModel.algKDTree, VicsekModel.VicsekModelVariables{Tf}, VicsekModel.VicsekModelParameters{Tf, Ti}}} where {Tf, Ti}","page":"Home","title":"VicsekModel.update!","text":"Update the position and velocity of particles using KDTree\n\nArguments\n\nvar: Variables of the Vicsek model\np: Parameters of the Vicsek model\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.update_position!-Union{Tuple{Tf}, Tuple{Matrix{Tf}, Matrix{Tf}, Vector{Tf}, Vector{Complex}, Tf, Tf, Base.OneTo}} where Tf","page":"Home","title":"VicsekModel.update_position!","text":"Update the position and velocity of particles\n\nArguments\n\npos: Position of particles\nvel: Velocity of particles\nθ: Angle of particles\nS: Complex number to calculate the average angle of neighbors\nv₀: Initial velocity\nη: Noise strength\nix: Index of particles\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.update_with_anim-Union{Tuple{Ti}, Tuple{Tf}, Tuple{Ti, VicsekModel.VicsekModelVariables{Tf}, VicsekModel.VicsekModelParameters{Tf, Ti}}} where {Tf, Ti}","page":"Home","title":"VicsekModel.update_with_anim","text":"Update the position and velocity of particles with animation\n\nArguments\n\nt: Time step\nvar: Variables of the Vicsek model\np: Parameters of the Vicsek model\n\n\n\n\n\n","category":"method"},{"location":"#VicsekModel.vicsekmodel_with_anim-Tuple{Any}","page":"Home","title":"VicsekModel.vicsekmodel_with_anim","text":"Vicsek model with animation\n\nArguments\n\nTmax: Maximum time step\n\n\n\n\n\n","category":"method"}]
}
