using PrecompileTools: @compile_workload

function workload()
    dt = DateTime(2021, 3, 28)
    igrf_Bd(1, 30, 4, dt)
    return
end

@compile_workload begin
    workload()
end
