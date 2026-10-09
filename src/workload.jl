using PrecompileTools: @compile_workload

function workload()
    dt = DateTime(2021, 3, 28)
    r = SM(1.0, 2.0, 3.0)
    sm2mag(r, dt)
    A = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    ts = [dt, dt + Day(1)]
    geo2gsm(A, ts)
    IGRF()(A, ts; from = GSM())
    get_mlt(A, ts)
    return
end

@compile_workload begin
    workload()
end
