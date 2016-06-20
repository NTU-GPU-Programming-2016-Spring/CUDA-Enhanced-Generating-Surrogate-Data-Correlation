
/* =============================================
                    Packages
   ============================================= */

// Native libraries below.
var fs         = require('fs');
// 3rd-party libraries below.
var app        = require('koa')();
var router     = require('koa-router')();
var static     = require('koa-static');
var bodyParser = require('koa-body');

/* =============================================
                 Web application
   ============================================= */

// Import the source of web application.
var web = require('./src/web.js');

// RESTful web path.
router
	.get ('/',      web.core.index)
	.post('/chart', web.chart.update);

// Enable the koa web server.
app
	// Set the static files path.
	.use(static('./static'))
	// Set the body parser.
	.use(bodyParser())
	// Set routes.
	.use(router.routes())
	.use(router.allowedMethods())
	// Start listening.
	.listen(8080);

// Error handling.
app
	.on('error', function (err) {
		log.error('server error', err);
	});