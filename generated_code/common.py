def generate_kernel_name_prefix(target):
  return f'{target}_' if target == 'gpu' else ''
