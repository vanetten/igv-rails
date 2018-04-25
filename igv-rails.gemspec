# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'igv/rails/version'

Gem::Specification.new do |spec|
  spec.name          = "igv-rails"
  spec.version       = Igv::Rails::VERSION
  spec.authors       = ["William Van Etten, PhD"]
  spec.email         = ["vanetten@bioteam.net"]
  spec.summary       = "IGV for Rails."
  spec.description   = "This gem provides igv javascript, css, and images for your Rails application."
  spec.homepage      = "https://github.com/vanetten/igv-rails"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0") + ["LICENSE", "README.md"]
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_dependency "railties", "~> 5.1"
  spec.add_development_dependency "bundler", "~> 1.16"
  spec.add_development_dependency "rake", "~> 12.3"
end
