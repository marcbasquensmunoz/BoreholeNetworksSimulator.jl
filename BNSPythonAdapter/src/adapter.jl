
BoreholeNetworksSimulator.unwrap(x::Py) = PythonCall.pyconvert(BoreholeNetworksSimulator.BoreholeOperation, x)

BoreholeNetworksSimulator.operate(operator::Py, i, options, Tfin, Tfout, Tb, q) = operator.operate(i, options, Tfin, Tfout, Tb, q)