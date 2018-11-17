"""
Change configuration file

Example:
    igdiscover config --set iterations 1

"""
import os
import logging
import ruamel.yaml

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
	pass


class Config:
	DEFAULT_PATH = 'igdiscover.yaml'

	def __init__(self, file):
		# Set some defaults.
		self.debug = False
		self.species = None
		self.sequence_type = 'Ig'
		self.merge_program = 'pear'
		self.flash_maximum_overlap = 300
		self.limit = None  # or an integer
		self.multialign_program = 'muscle-fast'
		self.minimum_merged_read_length = 300
		self.mismatch_penalty = None
		self.barcode_length = 0
		self.barcode_consensus = True
		self.iterations = 1
		self.ignore_j = False
		self.d_coverage = 70
		self.subsample = 1000
		self.stranded = False
		self.forward_primers = None
		self.reverse_primers = None
		self.rename = False
		self.race_g = False
		self.seed = 1
		self.exact_copies = None
		self.preprocessing_filter = dict(v_coverage=90, j_coverage=60, v_evalue=1E-3)
		self.pre_germline_filter = dict(
			unique_cdr3s=2,
			unique_js=2,
			whitelist=True,
			cluster_size=0,
			allow_stop=True,
			# allow_chimeras=False,
			cross_mapping_ratio=0.02,
			clonotype_ratio=0.12,
			exact_ratio=0.12,
			cdr3_shared_ratio=0.8,
			unique_d_ratio=0.3,
			unique_d_threshold=10
		)
		self.germline_filter = dict(
			unique_cdr3s=5,
			unique_js=3,
			whitelist=True,
			cluster_size=100,
			allow_stop=False,
			# allow_chimeras=False,
			cross_mapping_ratio=0.02,
			clonotype_ratio=0.12,
			exact_ratio=0.12,
			cdr3_shared_ratio=0.8,
			unique_d_ratio=0.3,
			unique_d_threshold=10
		)
		self.j_discovery = dict(allele_ratio=0.2, cross_mapping_ratio=0.1, propagate=True)
		self.cdr3_location = 'detect'

		self.read_from(file)

	def read_from(self, file):
		content = file.read()
		new_config = self.make_compatible(ruamel.yaml.safe_load(content))
		for key in ('preprocessing_filter', 'pre_germline_filter', 'germline_filter', 'j_discovery'):
			if key in new_config:
				self.__dict__[key].update(new_config[key])
				del new_config[key]
		self.__dict__.update(new_config)

	@staticmethod
	def make_compatible(config):
		"""
		Convert old-style configuration to new style.
		Raise ConfigurationError if configuration is invalid.
		Return updated config dict.
		"""
		if 'barcode_length' in config:
			raise ConfigurationError(
				'Old-style configuration of barcode length via "barcode_length"'
				'is no longer supported.')
		barcode_length_5prime = config.get('barcode_length_5prime', 0)
		barcode_length_3prime = config.get('barcode_length_3prime', 0)
		if barcode_length_5prime > 0 and barcode_length_3prime > 0:
			raise ConfigurationError(
				'barcode_length_5prime and barcode_length_3prime can currently '
				'not both be greater than zero.')
		if barcode_length_5prime > 0:
			config['barcode_length'] = barcode_length_5prime
		elif barcode_length_3prime > 0:
			config['barcode_length'] = -barcode_length_3prime
		config.pop('barcode_length_5prime', None)
		config.pop('barcode_length_3prime', None)

		if 'seed' in config and config['seed'] is False:
			config['seed'] = None

		for key in ('germline_filter', 'pregermline_filter'):
			if key in config and 'allele_ratio' in config[key]:
				config[key]['clonotype_ratio'] = config[key]['allele_ratio']
		return config

	@classmethod
	def from_default_path(cls):
		with open(cls.DEFAULT_PATH) as f:
			return Config(file=f)


class GlobalConfig:
	def __init__(self):
		self.use_cache = False
		path = os.getenv('XDG_CONFIG_HOME', os.path.expanduser('~/.config'))
		path = os.path.join(path, 'igdiscover.conf')
		if os.path.exists(path):
			with open(path) as f:
				config = ruamel.yaml.safe_load(f)
			if config is None:
				return
			self.use_cache = config.get('use_cache', False)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--set', nargs=2, default=[], metavar=('KEY', 'VALUE'), action='append',
		help='Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys.')
	arg('--file', default=Config.DEFAULT_PATH,
		help='Configuration file to modify. Default: igdiscover.yaml in current directory.')


def main(args):
	if args.set:
		with open(args.file) as f:
			config = ruamel.yaml.load(f, ruamel.yaml.RoundTripLoader)
		for k, v in args.set:
			v = ruamel.yaml.safe_load(v)
			# config[k] = v
			item = config
			# allow nested keys
			keys = k.split('.')
			for i in keys[:-1]:
				item = item[i]
			item[keys[-1]] = v
		tmpfile = args.file + '.tmp'
		with open(tmpfile, 'w') as f:
			print(ruamel.yaml.dump(config, Dumper=ruamel.yaml.RoundTripDumper), end='', file=f)
		os.rename(tmpfile, args.file)
	else:
		with open(args.file) as f:
			config = ruamel.yaml.safe_load(f)
		print(ruamel.yaml.dump(config), end='')
