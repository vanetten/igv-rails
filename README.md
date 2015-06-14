# igv-rails

## Installation

Add this line to the application's Gemfile:

```ruby
gem 'igv-rails'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install igv-rails

## Usage

Add IGV (and dependencies) to the application's Gemfile:

```ruby
gem 'font-awesome-rails'
gem 'jquery-rails'
gem 'jquery-ui-rails'
gem 'igv-rails'
```

Require IGV (and dependencies) to the application's application.js:

```javascript
//= require jquery
//= require jquery-ui
//= require igv
```

Require IGV CSS (and dependencies) to the application's application.css:

```css
/*
*= require font-awesome
*= require jquery-ui
*= require igv
*/
```

## Example

Provide a div container within a view:

```html
<div id="myDiv"></div>
```

Provide javascript to configure and load IGV within the view:

```javascript
$(document).ready(function () {
var div = $("#myDiv")[0],
options = {
    showNavigation: true,
    genome: "hg19",
    locus: "chr1:155,172,193-155,172,564",
    tracks: [
        {
            url: '//www.broadinstitute.org/igvdata/1KG/b37/data/NA06984/alignment/NA06984.mapped.ILLUMINA.bwa.CEU.low_coverage.20120522.bam',
            label: 'NA06984'
        }
    ]
};

igv.createBrowser(div, options);
});
```

## Contributing

1. Fork it ( https://github.com/[my-github-username]/igv-rails/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request

## Build Notes

This is [igv-web](https://www.broadinstitute.org/software/igv/home) GEMified for the Rails >= 4.1 asset pipeline through the following:

		bundle gem igv-rails
		cd igv-rails
		mkdir -p vendor/assets/javascripts
		curl http://www.broadinstitute.org/igv/projects/igv-web/dist/igv-all.min.js -o vendor/assets/javascripts/igv-all.js
		mkdir -p vendor/assets/stylesheets
		curl http://www.broadinstitute.org/igv/projects/igv-web/css/igv.css -o vendor/assets/stylesheets/igv.css
		mkdir -p vendor/assets/images
		curl http://www.broadinstitute.org/igv/projects/igv-web/dist/igv-all.min.map -o vendor/assets/images/igv-all.min.map
		echo "" >> README.md; echo "# igv appended README #" >> README.md; echo "" >> README.md
		curl https://github.com/broadinstitute/igv-web/blob/master/README.md >> README.md
		echo "" >> LICENSE; echo "# igv appended LICENSE #" >> LICENSE; echo "" >> LICENSE
		curl https://github.com/broadinstitute/igv-web/blob/master/license.txt >> LICENSE
		git add .
		git commit -am "igv-rails"
		git remote add origin git@github.com:vanetten/igv-rails.git

* modify **lib/igv/rails/version.rb** to match igv-all.js version

		VERSION = "0.0.1.*"

* modify **lib/igv/rails.rb** to subclass Rails::Engine

		class Engine < ::Rails::Engine
		end

* modify **igv-rails.gemspec**

		spec.summary       = "IGV for Rails."
		spec.description   = "This gem provides igv javascript, css, and images for your Rails application."
		spec.homepage      = "https://github.com/vanetten/igv-rails"
		spec.files         = `git ls-files -z`.split("\x0") + ["LICENSE", "README.md"]
		spec.add_dependency "railties", "~> 4.1"

* build

		rake build

* release

		rake release

## igv appended README

igv-web
=======

Lightweight HTML-5 versison of the Integrative Genomics Viewer (http://www.broadinstitute.org/igv).