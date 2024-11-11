const Router = require('koa-router');
const users = require('./users');
const valorations = require('./valorations')
const ingredients = require('./ingredientes');
const authRoutes = require('./authentication');
const dotenv = require('dotenv');
const jwtMiddleware = require('koa-jwt')
const scopeProtectedRoutes = require('./scopeExample.js');
const requests = require('./requests')



dotenv.config();


const router = new Router();
router.use('/valorations', valorations.routes());
router.use('/ingredientes', ingredients.routes());
router.use('/requests', requests.routes());



router.use(authRoutes.routes());


// desde esta linea todas las rutas requeriran un JWT, esto no aplica para
// las lineas anteriores
router.use(jwtMiddleware( { secret: process.env.JWT_SECRET } ))
router.use('/users', users.routes());
router.use('/scope-example', scopeProtectedRoutes.routes())

module.exports = router;