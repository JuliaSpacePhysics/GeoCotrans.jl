using PrecompileTools: @compile_workload

function workload()
    dt = DateTime(2021, 3, 28)
    B_true = (-46077.31133522, -14227.12618499, 233.14355744)
    all(igrf_B(SPH(6500, 30, 4), dt) .≈ B_true)
end

@compile_workload begin
    workload()
end
