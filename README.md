# Igv::Rails

This is [igv-web](https://www.broadinstitute.org/software/igv/home) GEMified for the Rails >= 4.1 asset pipeline through the following:

		bundle gem igv-rails
		cd igv-rails
		mkdir -p vendor/assets/javascripts
		curl http://www.broadinstitute.org/igv/projects/igv-web/dist/igv-all.min.js -o vendor/assets/javascripts/igv-all.min.js
		mkdir -p vendor/assets/stylesheets
		curl http://www.broadinstitute.org/igv/projects/igv-web/css/igv.css -o vendor/assets/stylesheets/igv.css
		echo "" >> README.md; echo "# igv appended README #" >> README.md; echo "" >> README.md
		curl https://github.com/broadinstitute/igv-web/blob/master/README.md >> README.md
		echo "" >> LICENSE; echo "# igv appended LICENSE #" >> LICENSE; echo "" >> LICENSE
		curl https://github.com/broadinstitute/igv-web/blob/master/license.txt >> LICENSE
		git add .
		git commit -am "igv-rails"
		git remote add origin git@github.com:vanetten/igv-rails.git

* modify **lib/igv/rails/version.rb** to match igv-all.min.js version

		VERSION = "0.0.1"

* modify **lib/igv/rails.rb** to subclass Rails::Engine

		class Engine < ::Rails::Engine
		end

* modify **igv-rails.gemspec**

	  spec.summary       = "IGV for Rails."
	  spec.description   = "This gem provides igv-all.min.js and igv.css for your Rails application."
	  spec.homepage      = "https://github.com/vanetten/igv-rails"
	  spec.files         = `git ls-files -z`.split("\x0") + ["LICENSE", "README.md"]
	  spec.add_dependency "railties", "~> 4.1"

* build

		rake build

* release

		rake release

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'igv-rails'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install igv-rails

## Usage

TODO: Write usage instructions here

## Contributing

1. Fork it ( https://github.com/[my-github-username]/igv-rails/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request

# igv appended README #

igv-web
=======

Lightweight HTML-5 versison of the Integrative Genomics Viewer (http://www.broadinstitute.org/igv).