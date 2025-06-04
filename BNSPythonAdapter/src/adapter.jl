
BoreholeNetworksSimulator.unwrap(x::Py) = PythonCall.pyconvert(BoreholeNetworksSimulator.BoreholeOperation, x)

BoreholeNetworksSimulator.operate(operator::Py, i, options, X) = operator.operate(i, options, X)