def generate_kernename_prefix(platform):
  return f'{platform}_' if platform == 'gpu' else ''
