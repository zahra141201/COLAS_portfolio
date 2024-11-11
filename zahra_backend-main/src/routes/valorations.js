const Router = require('koa-router');
const { Valoration } = require('../models');
const { Request } = require('../models');
const router = new Router();





router.post('/', async (ctx) => {
  try {
    console.log('Solicitud recibida para la creación de una valoracion');
    console.log('Datos de la solicitud:', ctx.request.body);

    const valoration = await Valoration.create(ctx.request.body);
    console.log('valoration creado con éxito:', valoration);

    ctx.status = 201;
    ctx.body = valoration;
  } catch (error) {
    console.error('Error al crear la valoration:', error);
    ctx.status = 500;
    ctx.body = { error: 'Ocurrió un error al crear la valoration' };
  }
});



// GET para obtener valoraciones por email_user
router.get('/:email_user', async (ctx) => {
  const { email_user } = ctx.params;

  try {
    // Buscar valoraciones por email_user
    const valorations = await Valoration.findAll({
      where: {
        email_user: email_user
      }
    });

    if (valorations.length === 0) {
      ctx.status = 404; // No encontrado
      ctx.body = { error: 'No se encontraron valoraciones para este usuario' };
      return;
    }

    ctx.status = 200; // OK
    ctx.body = valorations;
  } catch (error) {
    console.error('Error al obtener valoraciones por email_user:', error);
    ctx.status = 500; // Error del servidor
    ctx.body = { error: 'Ocurrió un error al obtener las valoraciones por email_user' };
  }
});





router.patch('/:id', async (ctx) => {
  const { id } = ctx.params;
  console.log(ctx.request.body); // Utilisez ctx.request.body au lieu de req.body
  const { comment, puntuation, email_user, made_by } = ctx.request.body;

  try {
    // Trouver la valoration par son ID
    const valoration = await Valoration.findByPk(id);

    if (!valoration) {
      ctx.status = 404;
      ctx.body = { error: 'Valoration not found' };
      return;
    }

    // Mettre à jour les champs de la valoration
    if (comment !== undefined) valoration.comment = comment;
    if (puntuation !== undefined) valoration.puntuation = puntuation;
    if (email_user !== undefined) valoration.email_user = email_user;
    if (made_by !== undefined) valoration.made_by = made_by;

    // Sauvegarder les changements dans la base de données
    await valoration.save();

    ctx.body = { message: 'Valoration updated successfully', valoration };
  } catch (error) {
    console.error('Error updating valoration:', error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while updating the valoration' };
  }
});

// Route pour supprimer une requête selon l'utilisateur et l'ingrédient
router.delete('/user/:made_by', async (ctx) => {
  try {
    const { made_by } = ctx.params;
    const valoration = await Valoration.findOne({
      where: { made_by },
    });

    if (valoration) {
      await valoration.destroy();
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


router.get('/', async (ctx) => {
  const { email_user, made_by } = ctx.query;

  try {
    // Vérifier s'il existe une évaluation pour email_user et made_by
    const existingValoration = await Valoration.findOne({
      where: {
        email_user: email_user,
        made_by: made_by
      }
    });

    if (existingValoration) {
      console.log('valoration creado con éxito:', existingValoration);
      ctx.status = 200; // OK
      ctx.body = existingValoration;
    } else {
      ctx.status = 404; // Non trouvé
      ctx.body = { error: 'Aucune évaluation trouvée pour cet utilisateur et cet évaluateur' };
    }
  } catch (error) {
    console.error('Erreur lors de la vérification de l\'existence de l\'évaluation :', error);
    ctx.status = 500; // Erreur du serveur
    ctx.body = { error: 'Une erreur s\'est produite lors de la vérification de l\'existence de l\'évaluation' };
  }
});



module.exports = router;