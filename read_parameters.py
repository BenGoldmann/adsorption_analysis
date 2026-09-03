def read_parameters(params_path):
	"""
	Very small parser for lines like:
	KEY = VALUE
	Ignores blank lines and comments starting with #. Returns a dictionary of key-value pairs.
	"""
	params = {}

	with open(params_path, "r") as f:
		for raw in f:
			line = raw.strip()
			if not line or line.startswith("#"):
				continue
			if "=" not in line:
				continue

			key, value = line.split("=", 1)
			params[key.strip()] = value.strip().strip('"').strip("'")
	return params