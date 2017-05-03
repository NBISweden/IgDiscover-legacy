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
		self.merge_program = 'pear'
		self.flash_maximum_overlap = 300
		self.limit = None  # or an integer
		self.cluster_program = 'vsearch'
		self.multialign_program = 'muscle-fast'
		self.maximum_expected_errors = None  # or an integer
		self.minimum_merged_read_length = 300
		self.mismatch_penalty = None
		self.barcode_length = 0
		self.iterations = 3
		self.ignore_j = False
		self.subsample = 1000
		self.stranded = False
		self.forward_primers = None
		self.reverse_primers = None
		self.rename = True
		self.race_g = False
		self.seed = 1
		self.exact_copies = None
		self.preprocessing_filter = dict(v_coverage=90, j_coverage=60, v_evalue=1E-3)
		self.pre_germline_filter = dict(unique_cdr3s=2, unique_js=2, check_motifs=False,
			whitelist=True, cluster_size=0, differences=0, allow_stop=True, cross_mapping_ratio=0.02,
			allele_ratio=0.1)
		self.germline_filter = dict(unique_cdr3s=5, unique_js=3, check_motifs=False,
			whitelist=True, cluster_size=100, differences=0, allow_stop=False, cross_mapping_ratio=0.02,
			allele_ratio=0.1)
		self.cdr3_location = [-80, -60]
		self.library_name = os.path.basename(os.getcwd())

		self.read_from(file)

	def read_from(self, file):
		content = file.read()
		new_config = self.make_compatible(ruamel.yaml.safe_load(content))
		for key in ('preprocessing_filter', 'pre_germline_filter', 'germline_filter'):
			if key in new_config:
				self.__dict__[key].update(new_config[key])
				del new_config[key]
		self.__dict__.update(new_config)

	def make_compatible(self, config):
		"""
		Convert old-style configuration to new style. Raise ConfigurationError if configuration is invalid.
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
		return config

	@classmethod
	def from_default_path(cls):
		with open(cls.DEFAULT_PATH) as f:
			return Config(file=f)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--set', nargs=2, default=[], metavar=('KEY', 'VALUE'), action='append',
		help='Set KEY to VALUE.')
	arg('--file', default=Config.DEFAULT_PATH,
		help='Configuration file to modify. Default: igdiscover.yaml in current directory.')


def main(args):
	if args.set:
		with open(args.file) as f:
			config = ruamel.yaml.load(f, ruamel.yaml.RoundTripLoader)
		for k, v in args.set:
			v = ruamel.yaml.safe_load(v)
			config[k] = v
		tmpfile = args.file + '.tmp'
		with open(tmpfile, 'w') as f:
			print(ruamel.yaml.dump(config, Dumper=ruamel.yaml.RoundTripDumper), end='', file=f)
		os.rename(args.file + '.tmp', args.file)
	else:
		with open(args.file) as f:
			config = ruamel.yaml.safe_load(f)
		print(ruamel.yaml.dump(config), end='')
