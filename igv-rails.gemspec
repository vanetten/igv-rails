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
  spec.description   = "This gem provides igv-all.min.js and igv.css for your Rails application."
  spec.homepage      = "https://github.com/vanetten/igv-rails"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0") + ["LICENSE", "README.md"]
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_dependency "railties", "~> 4.1"
  spec.add_dependency "font-awesome-rails", "~> 4.3"
  spec.add_dependency "jquery-rails", "~> 3.1"
  spec.add_dependency "jquery-ui-rails", "~> 5.0"
  spec.add_development_dependency "bundler", "~> 1.7"
  spec.add_development_dependency "rake", "~> 10.0"
end
