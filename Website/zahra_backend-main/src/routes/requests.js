const Router = require('koa-router');
const router = new Router();
const { Request } = require('../models');

// Route pour créer une requête
router.post('/', async (ctx) => {
  try {
    const { pick_up_date, comment, state, made_by, id_ingrediente } = ctx.request.body;
    const newRequest = await Request.create({
      pick_up_date,
      comment,
      state,
      made_by,
      id_ingrediente,
    });
    ctx.status = 201;
    ctx.body = newRequest;
  } catch (error) {
    console.error('Error creating request:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while creating the request' };
  }
});

// Route pour obtenir les requêtes selon un ingrédient
router.get('/ingredient/:id_ingrediente', async (ctx) => {
  try {
    const { id_ingrediente } = ctx.params;
    const requests = await Request.findAll({
      where: { id_ingrediente },
    });
    ctx.status = 200;
    ctx.body = requests;
  } catch (error) {
    console.error('Error fetching requests by ingredient:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while fetching the requests' };
  }
});

router.patch('/ingredient/:id_ingrediente', async (ctx) => {
  try {
    const { id_ingrediente } = ctx.params;
    const { id, pick_up_date, comment, state, made_by } = ctx.request.body;

    const request = await Request.findOne({
      where: { id_ingrediente, id }
    });

    if (request) {
      request.pick_up_date = pick_up_date;
      request.comment = comment;
      request.state = state;
      request.made_by = made_by;
      await request.save();
      ctx.status = 200;
      ctx.body = request;
    } else {
      ctx.status = 404;
      ctx.body = { error: 'Request not found' };
    }
  } catch (error) {
    console.error('Error updating request status:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while updating the request status' };
  }
});

// Route pour obtenir les requêtes selon un utilisateur
router.get('/user/:made_by', async (ctx) => {
  try {
    const { made_by } = ctx.params;
    const requests = await Request.findAll({
      where: { made_by },
    });
    ctx.status = 200;
    ctx.body = requests;
  } catch (error) {
    console.error('Error fetching requests by user:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while fetching the requests' };
  }
});

// Route pour supprimer une requête selon l'utilisateur et l'ingrédient
router.delete('/user/:made_by/ingredient/:id_ingrediente', async (ctx) => {
  try {
    const { made_by, id_ingrediente } = ctx.params;
    const request = await Request.findOne({
      where: { made_by, id_ingrediente },
    });

    if (request) {
      await request.destroy();
      ctx.status = 204; // No Content
    } else {
      ctx.status = 404;
      ctx.body = { error: 'Request not found' };
    }
  } catch (error) {
    console.error('Error deleting request:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while deleting the request' };
  }
});

module.exports = router;
