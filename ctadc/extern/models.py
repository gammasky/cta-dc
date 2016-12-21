"""
Model classes to generate XML.

TODO (about a week of work):
- add SourceLibrary class
- implement SourceLibrary.parse_xml()
- integrate this the existing Gammapy model classes to make analysis possible.
- don't couple this to gamma-cat. Gamma-cat should change to be natively in this format.
- sub-class Astropy Parameter and ParameterSet classes instead of starting from scratch?
- implement spatial and spectral mode registries instead of `if-elif` set on type to make SourceLibrary extensible.
- write test and docs
- Once modeling setup OK, ask new people to add missing models (see Gammalib, Fermi ST, naima, Sherpa, HESS)
  (it's one of the simplest and nicest things to get started with.

For XML model format definitions, see here:

* http://cta.irap.omp.eu/ctools/user_manual/getting_started/models.html#spectral-model-components
* http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/source_models.html
"""

__all__ = [
    'Parameter',
    'ParameterSet',
    'SourceModel',

    'SpectralModel',
    'SpectralModelPowerLaw',
    'SpectralModelPowerLaw2',
    'SpectralModelExpCutoff',

    'SpatialModel',
    'SpatialModelPoint',
    'SpatialModelGauss',
    'SpatialModelShell',

    'UnknownModelError',
    'MissingGammaCatDataError',
]


class UnknownModelError(ValueError):
    """
    Error when encountering unknown model types.
    """


class MissingGammaCatDataError(ValueError):
    """Exception to handle missing data in gamma-cat."""


class Parameter:
    def __init__(self, name, value, unit=''):
        self.name = name
        self.value = value
        self.unit = unit

    @classmethod
    def from_dict(cls, data):
        return cls(
            name=data['name'],
            value=data['val'],
            unit=data['unit'],
        )

    def to_xml(self):
        return '        <parameter name="{name}" value="{value}" unit="{unit}"/>'.format(**self.__dict__)


class ParameterSet:
    def __init__(self, data):
        """
        data : list of Parameter
        """
        self.data = data

    @classmethod
    def from_list_of_dict(cls, data):
        return cls([Parameter.from_dict(_) for _ in data])

    def par(self, name):
        """Access parameter by name"""
        for par in self.data:
            if name == par.name:
                return par

        raise IndexError('Parameter {} not found for pset: {}'.format(name, self.data))


class SourceModel:
    def __init__(self, source_name, source_type,
                 spatial_model, spectral_model):
        self.source_name = source_name
        self.source_type = source_type
        self.spatial_model = spatial_model
        self.spectral_model = spectral_model

    @classmethod
    def from_gammacat(cls, source):
        source_name = source.data['common_name']
        spectral_model = SpectralModel.from_gammacat(source)
        spatial_model = SpatialModel.from_gammacat(source)
        source_type = spatial_model.xml_source_type
        return cls(
            source_name=source_name,
            source_type=source_type,
            spatial_model=spatial_model,
            spectral_model=spectral_model,
        )

    def to_xml(self):
        data = dict(
            source_name=self.source_name,
            source_type=self.source_type,
            spectrum_xml=self.spectral_model.to_xml(),
            spatial_xml=self.spatial_model.to_xml(),
        )
        return """\
<source name="{source_name}" type="{source_type}">
    {spectrum_xml}
    {spatial_xml}
</source>
""".format(**data)


class SpectralModel:
    def __init__(self, parameters):
        self.parameters = parameters

    @property
    def parameter_xml(self):
        xml = [_.to_xml() for _ in self.parameters]
        return '\n'.join(xml)

    @classmethod
    def from_gammacat(cls, source):
        try:
            data = source.spectral_model.to_dict()
        except ValueError:
            raise MissingGammaCatDataError(source)

        pset = ParameterSet.from_list_of_dict(data['parameters'])
        if data['name'] == 'PowerLaw':
            model = SpectralModelPowerLaw.from_pset(pset)
        elif data['name'] == 'PowerLaw2':
            model = SpectralModelPowerLaw2.from_pset(pset)
        elif data['name'] == 'ExponentialCutoffPowerLaw':
            model = SpectralModelExpCutoff.from_pset(pset)
        else:
            print(pset)
            raise UnknownModelError('Unknown spectral model: {}'.format(data))

        return model

    def to_xml(self):
        return '<spectrum type="{}">\n{}\n    </spectrum>'.format(self.xml_name, self.parameter_xml)


class SpectralModelPowerLaw(SpectralModel):
    xml_name = 'PowerLaw'

    @classmethod
    def from_pset(cls, pset):
        par = pset.par('amplitude')
        prefactor = Parameter(name='Prefactor', value=par.value, unit=par.unit)
        par = pset.par('index')
        index = Parameter(name='Index', value=-par.value, unit=par.unit)
        par = pset.par('reference')
        scale = Parameter(name='Scale', value=par.value, unit=par.unit)

        parameters = [prefactor, index, scale]
        return cls(parameters=parameters)


class SpectralModelPowerLaw2(SpectralModel):
    xml_name = 'PowerLaw2'

    @classmethod
    def from_pset(cls, pset):
        par = pset.par('amplitude')
        integral = Parameter(name='Integral', value=par.value, unit=par.unit)
        par = pset.par('index')
        index = Parameter(name='Index', value=-par.value, unit=par.unit)
        par = pset.par('emin')
        lower_limit = Parameter(name='LowerLimit', value=par.value, unit=par.unit)
        par = pset.par('emax')
        upper_limit = Parameter(name='UpperLimit', value=par.value, unit=par.unit)

        parameters = [integral, index, lower_limit, upper_limit]
        return cls(parameters=parameters)


class SpectralModelExpCutoff(SpectralModel):
    xml_name = 'ExpCutoff'

    @classmethod
    def from_pset(cls, pset):
        par = pset.par('amplitude')
        prefactor = Parameter(name='Prefactor', value=par.value, unit=par.unit)
        par = pset.par('index')
        index = Parameter(name='Index', value=par.value, unit=par.unit)
        par = pset.par('lambda_')
        cutoff = Parameter(name='Cutoff', value=1 / par.value, unit=par.unit)
        par = pset.par('reference')
        scale = Parameter(name='Scale', value=par.value, unit=par.unit)

        parameters = [prefactor, index, cutoff, scale]
        return cls(parameters=parameters)


class SpatialModel:
    def __init__(self, parameters):
        self.parameters = parameters

    @property
    def parameter_xml(self):
        xml = [_.to_xml() for _ in self.parameters]
        return '\n'.join(xml)

    @classmethod
    def from_gammacat(cls, source):
        data = source.data

        if data['morph_type'] == 'point':
            model = SpatialModelPoint.from_gammacat(source)
        elif data['morph_type'] == 'gauss':
            model = SpatialModelGauss.from_gammacat(source)
        elif data['morph_type'] == 'shell':
            model = SpatialModelShell.from_gammacat(source)
        else:
            print(source.data)
            raise UnknownModelError('Unknown spatial model: {}'.format(source))

        return model

    def to_xml(self):
        return '<spatialModel type="{}">\n{}\n    </spatialModel>'.format(self.xml_spatial_type, self.parameter_xml)


class SpatialModelPoint(SpatialModel):
    xml_source_type = 'PointSource'
    xml_spatial_type = 'Point'

    @classmethod
    def from_gammacat(cls, source):
        d = source.data
        glon = Parameter(name='GLON', value=d['glon'].value, unit=d['glon'].unit)
        glat = Parameter(name='GLAT', value=d['glat'].value, unit=d['glat'].unit)

        parameters = [glon, glat]
        return cls(parameters=parameters)


class SpatialModelGauss(SpatialModel):
    xml_source_type = 'ExtendedSource'
    xml_spatial_type = 'Gauss'

    @classmethod
    def from_gammacat(cls, source):
        d = source.data
        glon = Parameter(name='GLON', value=d['glon'].value, unit=d['glon'].unit)
        glat = Parameter(name='GLAT', value=d['glat'].value, unit=d['glat'].unit)
        sigma = Parameter(name='Sigma', value=d['morph_sigma'].value, unit=d['morph_sigma'].unit)

        # TODO: fill `morph_sigma2` and `morph_pa` info

        parameters = [glon, glat, sigma]
        return cls(parameters=parameters)


class SpatialModelShell(SpatialModel):
    xml_source_type = 'ExtendedSource'
    xml_spatial_type = 'Shell'

    @classmethod
    def from_gammacat(cls, source):
        d = source.data
        glon = Parameter(name='GLON', value=d['glon'].value, unit=d['glon'].unit)
        glat = Parameter(name='GLAT', value=d['glat'].value, unit=d['glat'].unit)
        radius = Parameter(name='Radius', value=d['morph_sigma'].value, unit=d['morph_sigma'].unit)
        width = Parameter(name='Width', value=0, unit='deg')

        # raise NotImplementedError('SHELL')
        # import IPython; IPython.embed(); 1 / 0

        parameters = [glon, glat, radius, width]
        return cls(parameters=parameters)
