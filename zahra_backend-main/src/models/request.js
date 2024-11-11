'use strict';
const { Model } = require('sequelize');

module.exports = (sequelize, DataTypes) => {
  class Request extends Model {
    static associate(models) {
      this.belongsTo(models.User, {
        foreignKey: 'made_by',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les demandes associées lorsque l'utilisateur est supprimé
      });
      this.belongsTo(models.Ingredient, {
        foreignKey: 'id_ingrediente',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les demandes associées lorsque l'ingrédient est supprimé
      });
    }
  }

  Request.init({
    pick_up_date: DataTypes.DATE,
    comment: DataTypes.STRING,
    state: DataTypes.STRING,
    made_by: DataTypes.STRING,
    id_ingrediente: DataTypes.INTEGER
  }, {
    sequelize,
    modelName: 'Request',
  });

  return Request;
};
