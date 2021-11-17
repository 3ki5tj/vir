import os, sys
import webapp2, jinja2

JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))


class MainPage(webapp2.RequestHandler):
  def get(self, fnhtml="index.html"):
    self.response.headers['Content-Type'] = 'text/html'

    template = JINJA_ENVIRONMENT.get_template(fnhtml)

    # the template values are no longer used
    # the code is kept here simply for references
    template_values = {
        'last_update':
        '<p>Last updated on Mar. 12, 2014',
    }
    self.response.write(template.render(template_values))


class FirstPage(MainPage):
  def get(self):
    MainPage.get(self, "first.html")


class AboutPage(webapp2.RequestHandler):
  def get(self):
    self.response.headers['Content-Type'] = 'text/plain'
    self.response.write('test page')


app = webapp2.WSGIApplication(
  [('/', MainPage),
   ('/index', MainPage), # aliases
   ('/first', FirstPage),
   ('/about', AboutPage)],
  debug=True)

