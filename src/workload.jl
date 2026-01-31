using PrecompileTools: @compile_workload

function workload()
    dt = DateTime(2021, 3, 28)
    r = SM(1.0, 2.0, 3.0, dt)
    igrf_Bd(1, 30, 4, dt)
    MAG(r)
    return
end

@compile_workload begin
    workload()
end
