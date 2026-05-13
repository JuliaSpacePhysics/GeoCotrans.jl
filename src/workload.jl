using PrecompileTools: @compile_workload

function workload()
    dt = DateTime(2021, 3, 28)
    r = SM(1.0, 2.0, 3.0)
    igrf_Bd(1, 30, 4, dt)
    sm2mag(r, dt)
    return
end

@compile_workload begin
    workload()
end
