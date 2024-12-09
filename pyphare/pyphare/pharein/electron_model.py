from . import global_vars


class IsothermalClosure(object):
    closure_name = "isothermal"

    def __init__(self, **kwargs):
        self.Te = kwargs.get("Te", IsothermalClosure._defaultTe())

    @staticmethod
    def _defaultTe():
        return 0.1

    def dict_path(self):
        return {"name/": IsothermalClosure.closure_name, "Te": self.Te}

    @staticmethod
    def name():
        return IsothermalClosure.closure_name



class PolytropicClosure(object):
    closure_name = "polytropic"

    def __init__(self, **kwargs):
        self.Te = kwargs.get("Te", PolytropicClosure._defaultTe())
        self.Gamma = kwargs.get("Gamma", PolytropicClosure._defaultGamma())

    @staticmethod
    def _defaultTe():
        return 0.1

    @staticmethod
    def _defaultGamma():
        return 1.66

    def dict_path(self):
        return {"name/": PolytropicClosure.closure_name, "Te": self.Te, "Gamma": self.Gamma}

    @staticmethod
    def name():
        return PolytropicClosure.closure_name



class ElectronModel(object):
    def __init__(self, **kwargs):
        if kwargs["closure"] == "isothermal":
            self.closure = IsothermalClosure(**kwargs)
        elif kwargs["closure"] == "polytropic":
            self.closure = PolytropicClosure(**kwargs)
        else:
            self.closure = None

        global_vars.sim.set_electrons(self)

    def dict_path(self):
        return [
            ("electrons/pressure_closure/" + k, v)
            for k, v in self.closure.dict_path().items()
        ]
