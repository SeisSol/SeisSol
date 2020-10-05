def generate_kernename_prefix(target):
  return f'{target}_' if target == 'gpu' else ''
